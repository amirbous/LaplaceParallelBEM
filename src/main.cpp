#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"
#include "../include/geo_math.hpp"


int main(int argc, char* argv[]) {
    using ValueType = float;
    std::string mesh_file = (argc > 1 ? argv[1] : "sphere");

    size_t n_vertices = 0, n_faces = 0;
    std::vector<Vertex<ValueType>> vertices;
    std::vector<Face<ValueType>> faces;



    float mesh_centroid[] = {0.0f, 0.0f, 0.0f};
    float cent[] = {0.0f, 0.0f, 0.0f};
    float area_i{0.0f}, R{0.0f};


    read_vtk<ValueType>(mesh_file, vertices, faces, n_vertices, n_faces);
    std::cout << "Mesh: " << mesh_file << ", " << n_vertices << " nodes, " << n_faces << " faces" << std::endl;

    // preallocate matrix and vectors
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n_faces, n_faces); // main matrix
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(n_faces); // looad vector
    Eigen::VectorXd q; // Solution vector


    // Boundary Conditions
    for (auto& v : vertices) { v.potential = 10;} 


    // compute the center of the object
    for (auto& v : vertices) {

        mesh_centroid[0] += v.x;
        mesh_centroid[1] += v.y;
        mesh_centroid[2] += v.z;
    }

    mesh_centroid[0] /= static_cast<float>(n_vertices);
    mesh_centroid[1] /= static_cast<float>(n_vertices);
    mesh_centroid[2] /= static_cast<float>(n_vertices);


    std::cout << "Started assembling the matrix" << std::endl;
    
    /*-----------------------------------------------------------------------------
                 main loop to assemble the matrix: TODO parallelize
    -----------------------------------------------------------------------------*/
    for (size_t i = 0; i < n_faces; ++i) {
        cent[0] = (faces[i].v1->x + faces[i].v2->x + faces[i].v3->x) / 3;
        cent[1] = (faces[i].v1->y + faces[i].v2->y + faces[i].v3->y) / 3;
        cent[2] = (faces[i].v1->z + faces[i].v2->z + faces[i].v3->z) / 3;

        for (size_t j = 0; j < n_faces; ++j) {
            if (i == j) {
                // compute area of face:
                G(i, i) = regularized_integral(faces[j].v1, faces[j].v2, faces[j].v3, mesh_centroid);

            } else {
                // evaluating the gaussian integral
                G(i, j) = gauss_integral(faces[j].v1, faces[j].v2, faces[j].v3, cent, mesh_centroid);
            }
        }
    }
    /*-----------------------------------------------------------------------------
    -----------------------------------------------------------------------------*/


    // constructing load vector ==> averaged over faces
    for (size_t i = 0; i < n_faces; ++i) {
        phi(i) = (faces[i].v1->potential + faces[i].v2->potential + faces[i].v3->potential) / 3.0;
    }
    std::cout << "finished assembling the matrix" << std::endl;

    // solving the system
    q = G.partialPivLu().solve(phi);

    std::cout << "Solved G*q = phi for charge distribution." << std::endl;

    // getting the actual densities: per node and not per face
    std::vector<int> node_counts(n_vertices, 0);
    for (auto& v : vertices) { v.density = 0; }

    for (size_t i = 0; i < n_faces; ++i) {
        auto& f = faces[i];
        node_counts[f.v1->id]++; node_counts[f.v2->id]++; node_counts[f.v3->id]++;
        f.v1->density += q(i); f.v2->density += q(i); f.v3->density += q(i);
    }
    for (size_t i = 0; i < n_vertices; ++i) {
        if (node_counts[i] > 0) { vertices[i].density /= node_counts[i]; }
    }

    write_vtu<ValueType>(mesh_file, vertices, faces, n_vertices, n_faces);
    return 0;
}