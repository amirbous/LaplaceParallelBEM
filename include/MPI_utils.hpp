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

#endif