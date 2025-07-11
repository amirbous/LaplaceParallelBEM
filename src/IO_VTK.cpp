#include <fstream>
#include <iostream>

#include "../include/geometry.hpp"

template
<typename T>
void read_vtk(const std::string problem_name, Vertex<T>** vertices, Face<T>** faces) {

	std::ifstream fstream;
	std::string fname = "geo/" + problem_name + ".vtk";
	
	fstream.open(fname);

	if (!fstream) {
		std::cerr << "Error: file '" << fname << "' could not be opened or does not exist." << std::endl;
		
	}

}

template void read_vtk<double>(const std::string problem_name, Vertex<double>**, Face<double>**);

