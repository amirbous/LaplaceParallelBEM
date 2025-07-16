#include <iostream>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"


#include <string>
#include <vector>

int main() {

	using ValueType = double;

	std::string problem_name{"sphere"};	

	std::cout << "Boundary Element Method Laplace" << std::endl;


	size_t nvertices{0};
	size_t nfaces{0};
	std::vector<Vertex<double>> vertices;
	std::vector<Face<double>> faces;

	
	read_vtk<ValueType>(problem_name, vertices, faces, nvertices, nfaces);


	std::cout << "Example vertices" << std::endl;
	for (int i = 0; i < 10; i++) {
		std::cout << "Vi_x: " << vertices[i].x << std::endl;
 	}

 	for (int i = 0; i < nvertices; i++) {
 		vertices[i].potential = vertices[i].x + vertices[i].y + vertices[i].z;
 	}

 	std::cout << "Number of faces: " << nfaces << std::endl;

 	write_vtu<ValueType>(problem_name, vertices, faces, nvertices, nfaces);


	return 0;

}
