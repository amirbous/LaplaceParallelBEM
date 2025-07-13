#include <fstream>
#include <iostream>

#include "../include/geometry.hpp"

template
<typename T>
void read_vtk(const std::string problem_name, Vertex<T>** vertices, Face<T>** faces) {

	size_t nvertices{0}, nfaces{0}, ncellsOffsets{0};

	float xi{0.0f}, yi{0.0f}, zi{0.0f};
	int v1_i{0}, v2_i{0}, v3_i{0}, face_i{0}, face_j{0};
	int offseti{0};


	std::ifstream fstream;
	std::string line_buffer, word_buffer;
	int val_buffer{0};
	std::string fname = "geo/" + problem_name + ".vtk";


	const uint8_t nHeaderLines = 4;


	fstream.open(fname);

	if (!fstream) {
		std::cerr << "Error: file '" << fname << "' could not be opened or does not exist." << std::endl;
		
	}
	
	// skip 3 first lines: useless info
	for (int i = 0; i < nHeaderLines; i++) {
		std::getline(fstream, line_buffer);
		std::cout << line_buffer << std::endl;

	}

	// skip first word
	fstream >> word_buffer;
	fstream >> nvertices;
	fstream >> word_buffer;

	*vertices = (Vertex<T>*) malloc(sizeof(Vertex<T>) * nvertices);



	for (size_t v_id = 0; v_id < nvertices; v_id++) {

		fstream >> xi;
		fstream >> yi;
		fstream >> zi;
		(*vertices)[v_id] = Vertex<T>(v_id, xi, yi, zi);

	}

	// parse and store cells offsets
	fstream >> word_buffer;

	fstream >> nfaces;
	fstream >> ncellsOffsets;

	*faces = (Face<T>*) malloc(sizeof(Face<T>) * nfaces);



	fstream >> word_buffer;
	fstream >> word_buffer;


	for (size_t offset_id = 0; offset_id < nfaces; offset_id++) {
		fstream >> offseti;

	}

	fstream >> word_buffer;
	fstream >> word_buffer;




	for (size_t face_id = 0; face_id < nfaces; face_id++ ) {


			fstream >> v1_i;
			fstream >> v2_i;
			fstream >> v3_i;

			(*faces)[face_id] = Face<T>((*vertices) + v1_i, 
										(*vertices) + v2_i, 
										(*vertices) + v3_i
								  		);

			fstream >> val_buffer;
			fstream >> val_buffer;
			fstream >> val_buffer;



	}



}



template void read_vtk<double>(const std::string problem_name, Vertex<double>**, Face<double>**);

