# LaplaceParallelBEM
Parallel Boundary Element Method implementation for 3D Laplace

#### Boudary element solver for electrostatic dielectric problems using the Boundary Element Method

- TODO: 
	- [ ] seperate builds for MPI and serial (OpenMP version)
	- [ ] Parallelize building matrix: OpenMP
    - [ ] MPI distributed matrix assembly

- Implemented
    - [x] Load and write geometry with densities correctly
    - [x] Ensure convergence
    - [x] Link ginkgo as solver

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
