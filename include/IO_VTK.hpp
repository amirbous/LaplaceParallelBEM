#ifndef IO_VTK_HPP
#define IO_VTK_HPP


#include <string>

template <typename T>
struct Vertex;

template <typename T>
struct Face;


template
<typename T>
void read_vtk(const std::string problem_name, Vertex<T>** vertices, Face<T>** faces);

#endif

