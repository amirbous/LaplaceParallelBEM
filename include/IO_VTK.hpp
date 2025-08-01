#ifndef IO_VTK_HPP
#define IO_VTK_HPP


#include <string>
#include <vector>

template <typename T>
struct Vertex;

template <typename T>
struct Face;


template
<typename T>
void read_vtk(const std::string problem_name, std::vector<Vertex<T>>& vertices,
				std::vector<Face<T>>& faces,
				size_t& nvertices, size_t& nfaces);
template
<typename T>
void write_vtu(const std::string problem_name, const std::vector<Vertex<T>>& vertices, const std::vector<Face<T>>& faces,
							size_t &nvertices, size_t &nfaces);

#endif

