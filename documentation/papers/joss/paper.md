---
title: 'Ocellaris: a DG FEM solver for free surface flows'
tags:
  - DG FEM
  - Navier-Stokes
  - free surface
  - multi-phase flow
  - incompressible
  - FEniCS
authors:
 - name: Tormod Landet
   orcid: 0000-0001-5070-8056
   affiliation: "1"
affiliations:
 - name: Department of Mathematics, University of Oslo
   index: 1
date: 6 February 2019
bibliography: paper.bib
---

# Summary

Free surface flows are found wherever two immiscible fluids come into contact, such as at the interface between water and air in the ocean. Simulations of high Reynolds number free surface flows are important for the design of coastal, bottom fixed, and floating structures exposed to ocean waves, as well as partially filled pipes and tanks. The high density difference between water and air poses problems for numerical approximation across the interface, and the highly non-linear behaviour of the free surface, which can break and overturn, makes separating computations into two different fluid domains difficult. As a model for free surface flows, Ocellaris solves the variable density incompressible Navier-Stokes equations with discontinuous density fields,

$$
\newcommand{\pdiff}[2]{\frac{\partial#1}{\partial#2}}
\newcommand{\T}[1]{{\boldsymbol{{#1}}}}
\newcommand{\vel}{\T{u}}
\newcommand{\grav}{\T{g}}
\begin{aligned}
\rho \left( \pdiff{\vel}{t} + (\vel\cdot\nabla) \vel \right) &= \nabla\cdot\mu\left(\nabla \vel + (\nabla\vel)^T\right) - \nabla p + \rho \grav,\\
\nabla\cdot \vel &= 0,\\
\pdiff{\rho}{t} + \vel\cdot\nabla \rho &= 0.
\end{aligned}
$$

Low-order finite volume methods (FVM) are currently the most popular methods for solving the above equations when simulating high Reynolds number free surface flows[^1]. Open source FVM codes include OpenFOAM [@weller_FOAM_1998] and Gerris [@popinet_gerris_2003]. Potential flow methods are also used for simulation of non-breaking ocean waves [@tong_numerical_2019], and ship-wave interaction [@kring_nonlinear_1997; @faltinsen_hispeed_2005], but these methods require that viscous effects and vorticity can be disregarded. Finite volume methods are able to include these effects, and can produce exactly incompressible velocity fields. By this we mean that the velocity is pointwise divergence-free, the velocity facet fluxes sum to zero for each cell, and the velocity facet flux is continuous between neighbouring cells. In the mentioned FVM programs, the free surface advection is implemented with the volume of fluid (VOF) method [@hirt_volume_1981], which requires that the advected fluid indicator function $c$ is bounded, $c\in[0,1]$. Numerical VOF advection methods rely upon divergence free velocity fields to ensure mass conservation and bounded transport.

[^1]: FINE/Marine, FLOW-3D, Fluent, Orca3D, SHIPFLOW XCHAP, and StarCCM+ are examples of proprietary FVM free surface flow solvers used in the industry.

Due to the piecewise constant discretisation, low-order FVM methods do not suffer from Gibbs oscillations near the large jump in density and momentum at the interface between the immiscible fluids, as long as an appropriate flux limiter is applied [@leonard_adjusted_1979]. But there are some downsides to using low-order methods. Increasing the approximation order is more computationally efficient than increasing the mesh density in areas where the solution is expected to be smooth, which for wave simulations is most of the domain away from the free surface. Implementing higher order FVM methods is complicated on unstructured meshes, due to the need for large reconstruction stencils in order to obtain higher order approximations. Higher order finite element methods (FEM), such as the discontinuous Galerkin method (DG FEM), uses higher order basis functions in each cell to overcome this problem. The discontinuous nature of DG FEM methods additionally allows more computation to be performed locally with less coupling of cells, which can be beneficial for the overall computational efficiency [see e.g. @Kubatko09; @Kirby12].

Ocellaris is an exactly incompressible Navier-Stokes solver for free surface flows with a DG FEM based numerical method that supports higher order finite elements and contains specially designed slope limiting stabilisation filters to be able to handle large density transitions [@velslopelim; @OcellarisDoc]. Ocellaris is implemented in Python and C++ with FEniCS [@FEniCSBook; @FEniCS1.5] as the backend for the mesh and finite element assembly. PETSc is used for solving the resulting linear systems [@petsc-user-ref; @petsc-efficient; @Dalcin2011; @davis2004algorithm; @hypre-web-page].

Ocellaris uses a YAML based input file format documented in the user guide available at [ocellaris.org](https://www.ocellaris.org). The mesh geometry can be defined directly in the input file for simple geometries, or it can be loaded from any file format supported by meshio [@meshio], as long as the unstructured meshes are simplicial; triangles in 2D and tetrahedra in 3D. Due to the flexible nature of the implementation, custom numerical models can be added to a simulation by referencing external Python modules from the input file. Ocellaris uses the XDMF file format [@xdmf] for visualisations and a custom HDF5 format for restart files [@hdf5].

![Cylinder in waves. Rendered in Blender after using Paraview to extract the free surface from an Ocellaris simulation.](https://www.ocellaris.org/figures/cylinder_in_waves.jpg)

# References

