#ifndef GEO_MATH_HPP
#define GEO_MATH_HPP

template <typename T>
struct Vertex;

template <typename T>
struct Face;



// Define PI for calculations
#ifndef PI
#define PI 3.14159265358979323846
#endif

// Gaussian quadrature points and weights for triangles (6 points)
// Sauter/Schwab coefficients
static const int nGauss = 6; // degrees of sampling over master triangle
static const double gw[6] = {0.223381589678011, 0.223381589678011, 0.223381589678011,
                             0.109951743655322, 0.109951743655322, 0.109951743655322};
static const double gp[6][2] = {{0.659027622374092, 0.231933368553031}, {0.659027622374092, 0.108488943134895},
                                {0.231933368553031, 0.659027622374092}, {0.231933368553031, 0.108488943134895},
                                {0.108488943134895, 0.659027622374092}, {0.108488943134895, 0.231933368553031}};


// compute the area of 3d triangle
template
<typename T>
float face_area(const Vertex<T>* v1, const Vertex<T>* v2, 
					const Vertex<T>* v3);

// gauss integral: actual integral for the evaluation of matrix coefficients
template
<typename T>
double gauss_integral(const Vertex<T>* v1,
                    const Vertex<T>* v2, const Vertex<T>* v3, float cent[3]);




// regularized integral in case the integral triangle collides with the same triangle
template
<typename T>
double regularized_integral(const Vertex<T>* v1,
                    const Vertex<T>* v2, const Vertex<T>* v3);
#endif

// tarik masdoud
// salwa fi maheb arrih
// 