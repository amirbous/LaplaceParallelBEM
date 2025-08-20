# LaplaceParallelBEM
Parallel Boundary Element Method implementation for 3D Laplace

#### Boudary element solver for electrostatic dielectric problems using the Boundary Element Method

- TODO: 
    - [ ] Complete distributed solve: not-feasable splitting of data not structured since the beginning
    - [ ] Ginkgo distributed solver 
- Features Implemented
    - [x] Parallelize building matrix: both OpenMP and MPI 
    - [x] Solve with Ginkgo OpenMP or CUDA executor
    - [x] Parallelize building matrix: OpenMP
    - [x] Parallelize building matrix: MPI 
    - [x] OpenMP only version solves completely 

#
#### Requirements
- C++ 17 standard (15 should also work, but build configuration set to 17)
- Ginkgo v1.10.0
- OpenMP
- MPI_CXX
- python v3.10+ only for automated tests



#### Build instructions 

- configure and build

Important: The user should provide a Ginkgo's installation directory `GINKGO_USR` when invoking the first cmake call 

example build command:

`cmake -S -DGINKGO_USR=<Gingko root directory> . -B build`

`cmake --build build -j`

Inside the build 2 main directores will form:

    - serial: has a serial version program `electromain`, which also uses OpenMP.

    -  parallel: has the compiled MPI program `electroparal`, can be run as a serial program without using MPI

- running the program

For both program, the geometry problem name `<prolem name>` should be specified as the first command line arguments. A corresonding file `<problem_name>`.vtk should exist. The current implementation can only work with VTK 5.0 2D second element order meshe files. The repo has already some examples of different models.

The serial version has 2 modes `<precompute centroids>` specified as 0 or 1 as the next command line argument after the geometry problem name (default is 0), which are either to precompute and store a certain auxiliary numerical attribute Centroids and access them from the array during the computation

To run each of the programs can be done using 

    - serial: `./electromain <problem name> <precompute centroids>` 

    - parallel: `mpirun -np <n_mpi_cores> ./electroparal <problem name>  

 
#### Testing
After building, it is possible to run `ctest` in the build directory, which will execute a series of tests

The automated tests don't cover the complete convergence of the method, as they are meant to be lightweight and check partial computations at correctness of the results at different stages of simulation.

### Convergence

The implementation has been compared against an external solver for converngece.

The model Crane joint has been used for this convergence test


### Results and scaling

- OpenMP

- MPI

- Hybird OpenMP-MPI

#### Example models

1. Cylinder:

<table>
  <tr>
    <td align="center">
      <b>Potential</b><br>
      <img src="https://github.com/amirbous/LaplaceParallelBEM/blob/96d5882a0502883bc281958e2b723c4858c93e5a/screenshots/CylinderPotential.png" alt="Potential" width="400"/>
    </td>
    <td align="center">
      <b>Electrostatic Field Density</b><br>
      <img src="https://github.com/amirbous/LaplaceParallelBEM/blob/96d5882a0502883bc281958e2b723c4858c93e5a/screenshots/CylinderDensity.png" alt="Density" width="400"/>
    </td>
  </tr>
</table>

2. Crane joint
<table>
  <tr>
    <td align="center">
      <b>Potential</b><br>
      <img src="https://github.com/amirbous/LaplaceParallelBEM/blob/b5b749fa03a002b703f8aec7aa5622f48a7d83ae/screenshots/TraxStickPotential.png" alt="Potential" width="400"/>
    </td>
    <td align="center">
      <b>Electrostatic Field Density</b><br>
      <img src="https://github.com/amirbous/LaplaceParallelBEM/blob/b5b749fa03a002b703f8aec7aa5622f48a7d83ae/screenshots/TraxStickDensity.png" alt="Density" width="400"/>
    </td>
  </tr>
</table>
