#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


template 
<typename T>
struct Vertex {
	Vertex (): x(0.0f), y(0.0f), z(0.0f), 
				potential(0.0), density(-1.0)  {}
	Vertex (float x, float y, float z): x(x), y(y), z(z),
				potential(0.0), density(-1.0){
	}
	Vertex (float x, float y, float z, T potential) : 
				x(x), y(y), z(z), 
				potential(potential), density(-1.0) {
			}
	Vertex (float x, float y, float z, T potential, T density) : 
				x(x), y(y), z(z), 
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
	Vertex<T> v1;
	Vertex<T> v2;
	Vertex<T> v3;

};

#endif
