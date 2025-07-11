#include <iostream>

#include "../include/geometry.hpp"
#include "../include/IO_VTK.hpp"

#include <string>

int main() {

	std::string problem_name{"sphere"};	

	std::cout << "Boundary Element Method Laplace" << std::endl;
	
	read_vtk<double>(problem_name, nullptr, nullptr);

	return 0;

}
