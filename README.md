# LaplaceParallelBEM
Parallel Boundary Element Method implementation for 3D Laplace

#### Boudary element solver for electrostatic dielectric problems using the Boundary Element Method

- TODO: 
    - [ ] Move Laplace-smoothing kernel to a seperate function (minor) 
    - [ ] Parallelize building matrix: OpenMP **In progress**
    - [ ] MPI distributed matrix assembly


- Implemented
    - [x] Load and write geometry with densities correctly
    - [x] Ensure convergence
    - [x] Laplace-smoothing kernel for more homogeneous less noisy solution profile
    - [x] Link ginkgo as solver

#### Example models
1. Cylinder:


