# Overview
Structural topology optimization is a mathematical methodology that helps in selecting a suitable initial shape for a new structure, facilitating and speeding up its design. The aim is to decide how the material should be distributed within a domain so that the structure has the greatest possible rigidity, supporting the application of external loads without suffering excessive displacements and deformations, remaining in static equilibrium and satisfying a maximum volume constraint.

This is a Matlab implementation of a robust and efficient algorithm for solving large-scale three-dimensional structural topology optimization problems, in which the optimization problem is solved by a globally convergent sequential linear programming (SLP) method with a stopping criterion based on first-order optimality conditions. The SLP approach is combined with a multiresolution scheme, that employs different discretizations to deal with displacement, design and density variables, allowing high-resolution structures to be obtained at a low computational cost. In addition, the multigrid method is employed as a preconditioner for the conjugate gradient method, substantially reducing the time spent solving the linear equilibrium systems. Since multiresolution can lead to the appearance of unwanted artefacts in the structure, it is possible to use an adaptive strategy for increasing the degree of approximation of the displacements, with a technique for suppressing unnecessary problem variables, which provides more accurate solutions without excessively reducing the efficiency of the algorithm. The implementation also has a thresholding strategy, based on gradient information, to obtain structures composed only by solid or void regions.

A scientific article containing some results obtained with this implementation is avaiable at [arxiv.org/abs/2311.05595](https://arxiv.org/abs/2311.05595).

# Guide
A guide on how to start using the program is detailed at the [User Manual](User-Manual-TopOptSLP3D.pdf). 

# License
This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

Credits must be given to the authors: Alfredo Vitorino and Francisco A. M. Gomes. 

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
