#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


struct Vertex {
	Vertex (): x(0.0f), y(0.0f), z(0.0f), 
	           potential(0.0f), density(-1.0f), id(0) {}

	Vertex (float x, float y, float z, int id): 
	       x(x), y(y), z(z), potential(0.0f), density(-1.0f), id(id) {}

	Vertex (float x, float y, float z, float potential, int id): 
	       x(x), y(y), z(z), potential(potential), density(-1.0f), id(id) {}

	Vertex (float x, float y, float z, float potential, float density, int id): 
	       x(x), y(y), z(z), potential(potential), density(density), id(id) {}

	Vertex(const Vertex& other): 
	       x(other.x), y(other.y), z(other.z),
	       potential(other.potential), density(other.density), id(other.id) {}

    bool operator<(const Vertex& other) const {
        return id < other.id;
    }


    bool operator==(const Vertex& other) const {
        return id == other.id;
    }

	float x;
	float y;
	float z;

	float potential;
	float density;

	int id;
};


struct Face {

	Face(int v1, int v2, int v3)
        : v1(v1), v2(v2), v3(v3) {}
    Face()
        : v1(0), v2(0), v3(0) {}

        
    int v1;
    int v2;
    int v3;



};


#endif
