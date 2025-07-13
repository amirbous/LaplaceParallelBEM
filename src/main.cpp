#include <iostream>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"

#include <string>

int main() {

	std::string problem_name{"sphere"};	

	std::cout << "Boundary Element Method Laplace" << std::endl;

	Vertex<double>* vertices = nullptr;
	Face<double>* faces = nullptr;
	
	read_vtk<double>(problem_name, &vertices, &faces);


	std::cout << "Example vertices" << std::endl;
	for (int i = 0; i < 10; i++) {
		std::cout << "Vi_x" << vertices[i].x << std::endl;
 	}

 	std::cout << "Example faces" << std::endl;
	for (int i = 0; i < 10; i++) {
		std::cout << "Fi_v1 " << faces[i].v1->id << std::endl;
 	}

	return 0;

}
