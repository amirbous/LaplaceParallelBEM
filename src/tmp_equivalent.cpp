#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <string>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"

constexpr double PI = 3.14159265358979323846;
static const int nGauss = 6;
static const double gw[6] = {0.223381589678011, 0.223381589678011, 0.223381589678011,
                             0.109951743655322, 0.109951743655322, 0.109951743655322};
static const double gp[6][2] = {{0.659027622374092, 0.231933368553031},
                                {0.659027622374092, 0.108488943134895},
                                {0.231933368553031, 0.659027622374092},
                                {0.231933368553031, 0.108488943134895},
                                {0.108488943134895, 0.659027622374092},
                                {0.108488943134895, 0.231933368553031}};

// Compute triangle normal (not normalized)
Eigen::Vector3d triangle_normal(const Vertex<double>& v1,
                                const Vertex<double>& v2,
                                const Vertex<double>& v3) {
    Eigen::Vector3d e1(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
    Eigen::Vector3d e2(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
    return e1.cross(e2);
}

// Compute centroid of triangle
Eigen::Vector3d triangle_centroid(const Vertex<double>& v1,
                                  const Vertex<double>& v2,
                                  const Vertex<double>& v3) {
    return Eigen::Vector3d(
        (v1.x + v2.x + v3.x) / 3.0,
        (v1.y + v2.y + v3.y) / 3.0,
        (v1.z + v2.z + v3.z) / 3.0
    );
}

// High-order quadrature over a triangle for kernels 1/r and dG/dn
void gauss_integral(const Eigen::Vector3d& xi,
                    const Eigen::Vector3d& n_src,
                    const Vertex<double>* v1,
                    const Vertex<double>* v2,
                    const Vertex<double>* v3,
                    double& Gij, double& Hij) {
    Eigen::Vector3d y1(v1->x, v1->y, v1->z);
    Eigen::Vector3d y2(v2->x, v2->y, v2->z);
    Eigen::Vector3d y3(v3->x, v3->y, v3->z);
    double area = 0.5 * triangle_normal(*v1, *v2, *v3).norm();
    Gij = 0.0; Hij = 0.0;
    for (int k=0; k<nGauss; ++k) {
        double u = gp[k][0], v = gp[k][1];
        if (u+v>1) { u = 1-u; v = 1-v; }
        Eigen::Vector3d y = y1 + u*(y2-y1) + v*(y3-y1);
        Eigen::Vector3d r = xi - y;
        double rnorm = r.norm();
        if (rnorm < 1e-15) continue;    
        double w = gw[k] * area;
        Gij += w * (1.0/(4.0*PI*rnorm));
        double dG = r.dot(n_src) / (4.0*PI*std::pow(rnorm,3));
        Hij += w * dG;
    }
}

int main(int argc, char* argv[]) {
    using ValueType = double;
    std::string problem = (argc>1 ? argv[1] : "sphere");

    size_t nV=0, nF=0;
    std::vector<Vertex<ValueType>> vertices;
    std::vector<Face<ValueType>> faces;
    read_vtk<ValueType>(problem, vertices, faces, nV, nF);
    std::cout << nV << " nodes, " << nF << " faces" << std::endl;

    // Dirichlet phi
    for (auto& v: vertices)
        //v.potential = v.x*v.x - v.y*v.y;
        v.potential = 0;

    size_t N = faces.size();
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N,N);
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(N);

    // Precompute centroids and normals
    std::vector<Eigen::Vector3d> cent(N), normals(N);
    std::vector<double> area(N);
    for (size_t i=0;i<N;++i) {
        cent[i] = triangle_centroid(*faces[i].v1,*faces[i].v2,*faces[i].v3);
        normals[i] = triangle_normal(*faces[i].v1,*faces[i].v2,*faces[i].v3);
        if (normals[i].dot(cent[i])<0) normals[i] = -normals[i];
        area[i] = 0.5 * normals[i].norm();
        normals[i].normalize();
    }

    // Assemble
    for (size_t i=0;i<N;++i) {
        phi(i) = (faces[i].v1->potential + faces[i].v2->potential + faces[i].v3->potential)/3.0;
        for (size_t j=0;j<N;++j) {
            if (i==j) {
                H(i,j) = 0.5;
                double Gii, Hii;
                gauss_integral(cent[i], normals[j], faces[j].v1, faces[j].v2, faces[j].v3, Gii, Hii);
                G(i,j) = Gii;
            } else {
                Eigen::Vector3d r = cent[i] - cent[j];
                double rnorm = r.norm();
                if (rnorm < 2*std::sqrt(area[j])) {
                    double Gij, Hij;
                    gauss_integral(cent[i], normals[j], faces[j].v1, faces[j].v2, faces[j].v3, Gij, Hij);
                    G(i,j)=Gij; H(i,j)=Hij;
                } else {
                    G(i,j) = (1.0/(4.0*PI*rnorm))*area[j];
                    H(i,j) = (r.dot(normals[j])/(4.0*PI*std::pow(rnorm,3))) * area[j];
                }
            }
        }
    }

    // Solve
    Eigen::VectorXd rhs = H * phi;
    Eigen::VectorXd q = G.colPivHouseholderQr().solve(rhs);

    // Error
    double sum2=0, A=0, mE=0;
    for (size_t i=0;i<N;++i) {
        Eigen::Vector3d grad(2.0*cent[i].x(), -2.0*cent[i].y(), 0.0);
        double qex = normals[i].dot(grad);
        double e = std::abs(q(i)-qex);
        sum2 += e*e * area[i]; A += area[i]; mE = std::max(mE,e);
    }
    double L2 = std::sqrt(sum2/A);
    std::cout<<"L2 error:"<<L2<<" max:"<<mE<<std::endl;

    // Visualize
    std::vector<int> cnt(nV,0);
    for (size_t i=0;i<N;++i) {
        auto& f = faces[i];
        cnt[f.v1->id]++; cnt[f.v2->id]++; cnt[f.v3->id]++;
        f.v1->density += q(i);
        f.v2->density += q(i);
        f.v3->density += q(i);
    }
    for (size_t i=0;i<nV;++i) if(cnt[i]) vertices[i].density /= cnt[i];

    write_vtu<ValueType>(problem, vertices, faces, nV, nF);
    return 0;
}
