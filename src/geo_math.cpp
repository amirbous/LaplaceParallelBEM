#include "../include/geometry.hpp"
#include "../include/geo_math.hpp"

#include <cmath> 


/*
*
*	This file includes all the geometrical math kernels needed
*
*/

// mesh centroid is also needed to ensure the correct sign of the normal 
template
<typename T>
float face_area(const Vertex<T>* v1, const Vertex<T>* v2, 
					const Vertex<T>* v3){
	
	// vector to store the normal vectir (x, y, z)
	float triangle_normal[3] = {0.0f};

	// unit vector of triangle 1 and 2 (e1, e2) and centroid
	float cent[3] = {0.0f};
	float e1[3] = {0.0f}, e2[3] = {0.0f}, f_normal[3] = {0.0f};

	// compute centroide
	cent[0] = (v1->x + v2->x + v3->x) / 3;
	cent[1] = (v1->y + v2->y + v3->y) / 3;
	cent[2] = (v1->z + v2->z + v3->z) / 3;


	// compute unit vectors in triangle coordinates
	e1[0] = v2->x - v1->x;
	e1[1] = v2->y - v1->y;
	e1[2] = v2->z - v1->z;


	e2[0] = v3->x - v1->x;
	e2[1] = v3->y - v1->y;
	e2[2] = v3->z - v1->z;

	// compute cross product and save in e1
	// normal is now stored in e1
	f_normal[0] = e1[1] * e2[2] - e1[2] * e2[1];
	f_normal[1] = e1[2] * e2[0] - e1[0] * e2[2];
	f_normal[2] = e1[0] * e2[1] - e1[1] * e2[0];



	float area = 0.5 * std::sqrt(f_normal[0] * f_normal[0]
							 + f_normal[1] * f_normal[1]
							 + f_normal[2] * f_normal[2]);

	return area;
}


template
<typename T>
double gauss_integral(const Vertex<T>* v1,
                    const Vertex<T>* v2, const Vertex<T>* v3,
                    float cent[3]) {
    

    double Gij = 0.0;
    // local coordinate vector variables
    double u{0.0}, v{0.0}, rnorm{0.0}, w{0.0};
    float area{0.0f};

    // vectors to handle all transformations
    float y1[] = {v1->x, v1->y, v1->z};
    float y2[] = {v2->x, v2->y, v2->z};
    float y3[] = {v3->x, v3->y, v3->z};

    // passage to double when evaluating 
    // coordinates of master triangle
    double y[] = {0.0, 0.0, 0.0};


    // variable to store the result
    area = face_area(v1, v2, v3);

    for (int k = 0; k < nGauss; k++) {

        // local coordinates in the ansatz space
        u = gp[k][0], v = gp[k][1];
        if (u + v > 1.0) { 
            u = 1.0 - u; v = 1.0 - v; 
        }

        // gaussian kernel for local coordinates
        // transitioning to double datatypes, accuracy imposed by gaussian kernel
        y[0] = y1[0] + u * (y2[0] - y1[0]) + v * (y3[0] - y1[0]) - cent[0];
        y[1] = y1[1] + u * (y2[1] - y1[1]) + v * (y3[1] - y1[1]) - cent[1];
        y[2] = y1[2] + u * (y2[2] - y1[2]) + v * (y3[2] - y1[2]) - cent[2];

        rnorm = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);

        if (rnorm < 1e-15) {
            continue;
        }

        w = gw[k] * area;
        Gij += w * (1.0 / (4.0 * PI * rnorm));

    }

    return Gij;

}

template
<typename T>
double regularized_integral(const Vertex<T>* v1,
                    const Vertex<T>* v2, const Vertex<T>* v3
                    ) {

    double Gii{0.0};
    float area{0.0f}, R{0.0f};

    area = face_area(v1, v2, v3);

    R = std::sqrt(area / PI);

    Gii = R / 2.0;

    return Gii;

}


template double gauss_integral<float>(const Vertex<float>* v1,
                    const Vertex<float>* v2, const Vertex<float>* v3,
                    float cent[3]);


template double gauss_integral<double>(const Vertex<double>* v1,
                    const Vertex<double>* v2, const Vertex<double>* v3,
                    float cent[3]);


template float face_area<float>(const Vertex<float>* v1, const Vertex<float>* v2, 
					const Vertex<float>* v3);
template float face_area<double>(const Vertex<double>* v1, const Vertex<double>* v2, 
					const Vertex<double>* v3);

template double regularized_integral<float>(const Vertex<float>* v1,
                    const Vertex<float>* v2, const Vertex<float>* v3);

template double regularized_integral<double>(const Vertex<double>* v1,
                    const Vertex<double>* v2, const Vertex<double>* v3);