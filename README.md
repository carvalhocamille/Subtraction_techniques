# Subtraction_techniques

In this repository you will find MATLAB codes that compute a modified representation of layer potentials to avoid the close evaluation problem.
Codes have been used to produce numerical results in the manuscript entitled "Modified representations for the close evaluation problem", Carvalho (2021).

* There are two folders: 2D for two dimensional problems, 3D for three dimensional problems
* In 2D:
  * Run `SingleLayerPotential_subtraction.m` to compute the modified representation of the solution of the exterior Neumann Laplace's problem
  * Run `Helmholtz_subtraction_techniques.m` to compute the modified representation of the solution of the sound-soft scattering problem
  * `BoundaryCurve.m` allows you to change the boundary 
* In 3D:
  * Run `SLP3Dbasic_density_subtraction.m` to compute the modified representation of the solution of the exterior Neumann Laplace's problem
  * Run `helmholtz3D.m` to compute the modified representation of the solution of the sound-soft scattering problem
  * `ComputeSurface.m` allows you to change the boundary 


### Funding Support
NSF Grant DMS-1819052 

### Licence MIT
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5523373.svg)](https://doi.org/10.5281/zenodo.5523373)


