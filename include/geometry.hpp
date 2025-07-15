#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


template 
<typename T>
struct Vertex {
	int id;

	Vertex (): id(0), x(0.0f), y(0.0f), z(0.0f), 
				potential(0.0), density(-1.0)  {}
	Vertex (int id, float x, float y, float z): 
				id(id), x(x), y(y), z(z),
				potential(0.0), density(-1.0){
	}
	Vertex (int id, float x, float y, float z, T potential) : 
				id(id), x(x), y(y), z(z), 
				potential(potential), density(-1.0) {
			}
	Vertex (int id, float x, float y, float z, T potential, T density) : 
				id(id), x(x), y(y), z(z), 
				potential(potential), density(density)  {
			}
	float x;
	float y;
	float z;

	T potential;
	T density;



};

template
<typename T>
struct Face {

	Face(Vertex<T>* v1, Vertex<T>* v2, Vertex<T>* v3)
        : v1(v1), v2(v2), v3(v3) {}
        
    Vertex<T>* v1;
    Vertex<T>* v2;
    Vertex<T>* v3;



};

#endif
