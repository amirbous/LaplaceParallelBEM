#ifndef IO_VTK_HPP
#define IO_VTK_HPP


#include <string>
#include <vector>

template <typename T>
struct Vertex;

struct Face;


template
<typename T>
void read_vtk(const std::string problem_name, std::vector<Vertex<T>>& vertices,
				std::vector<Face>& faces,
				int& nvertices, int& nfaces);
template
<typename T>
void write_vtu(const std::string problem_name, const std::vector<Vertex<T>>& vertices, const std::vector<Face>& faces,
							int &nvertices, int &nfaces);

#endif

