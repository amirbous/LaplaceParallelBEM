#ifndef IO_VTK_HPP
#define IO_VTK_HPP


#include <string>
#include <vector>

struct Vertex;

struct Face;



void read_vtk(const std::string problem_name, std::vector<Vertex>& vertices,
				std::vector<Face>& faces,
				int& nvertices, int& nfaces);

void write_vtu(const std::string problem_name, const std::vector<Vertex>& vertices, const std::vector<Face>& faces,
							int &nvertices, int &nfaces);

#endif

