# LaplaceParallelBEM
Parallel Boundary Element Method implementation for 3D Laplace

#### Boudary element solver for electrostatic dielectric problems using the Boundary Element Method

- TODO: 
    - [ ] Parllelize other kernels with OpenMP: density mappings and Laplce smoothing - Will not result in too much speed up but interesting parallelization problems (reduction over array)
    - [ ] Complete distributed solve: not-feasable splitting of data not structured since the beginning
    - [ ] Maybe: GInkgo distributed solver 
- Implemented
    - [x] Parallelize building matrix: both OpenMP and MPI 
    - [x] Solve with Ginkgo OpenMP or CUDA executor
    - [x] Parallelize building matrix: OpenMP
    - [x] Parallelize building matrix: MPI 
    - [x] OpenMP only version solves completely 
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
