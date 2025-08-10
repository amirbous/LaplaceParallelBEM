#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <mpi.h>

MPI_Datatype MPI_FACE;

void create_MPIFace_type() {
    Face face;

    int block_lengths[3] = {1, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};

    MPI_Aint base_address;
    MPI_Get_address(&face, &base_address);
    MPI_Get_address(&face.v1, &displacements[0]);
    MPI_Get_address(&face.v2, &displacements[1]);
    MPI_Get_address(&face.v3, &displacements[2]);

    for (int i = 0; i < 3; i++) {
        displacements[i] = displacements[i] - base_address;
    }

    MPI_Type_create_struct(3, block_lengths, displacements, types, &MPI_FACE);
    MPI_Type_commit(&MPI_FACE);
}


MPI_Datatype MPI_VERTEX;

void create_MPIVertex_type() {
    Vertex vertex;

    int block_lengths[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[6];
    MPI_Datatype types[6] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};

    MPI_Aint base_address;
    MPI_Get_address(&vertex, &base_address);
    MPI_Get_address(&vertex.x, &displacements[0]);
    MPI_Get_address(&vertex.y, &displacements[1]);
    MPI_Get_address(&vertex.z, &displacements[2]);
    MPI_Get_address(&vertex.potential, &displacements[3]);
    MPI_Get_address(&vertex.density, &displacements[4]);
    MPI_Get_address(&vertex.id, &displacements[5]);

    for (int i = 0; i < 6; i++) {
        displacements[i] -= base_address;
    }

    MPI_Type_create_struct(6, block_lengths, displacements, types, &MPI_VERTEX);
    MPI_Type_commit(&MPI_VERTEX);
}


#endif

