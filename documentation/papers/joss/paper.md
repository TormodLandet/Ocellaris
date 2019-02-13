---
title: 'Ocellaris: an exactly incompressible DG FEM solver for free surface flows'
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

Free surface flows are found wherever two immiscible fluids come into contact, such as at the interface between water and air in the ocean. Simulation of high Reynolds number free surface flows are important for the design of coastal, bottom fixed, and floating structures exposed to ocean waves, as well as partially filled pipes and tanks. The high density difference between water and air poses problems for numerical approximations of both fluids simultaneously, and the highly non-linear behaviour of free surface, which can break and overturn, makes separating computations into two different fluid domains difficult.

Low order finite volume methods are currently what is most used for solving the Navier-Stokes equations when simulating high Reynolds number free surface flows, though boundary element potential flow methods can also give reasonable results for non-breaking ocean waves and for simulation of ship-wave interaction [@weller_FOAM_1998; @popinet_gerris_2003; @kring_nonlinear_1997]. Finite volume methods can produce exactly incompressible velocity fields, which are important for the advection of the free surface, and, due to the piecewise constant discretisation, they do not suffer from Gibbs oscillations near the large jump in density and momentum at the interface between the immiscible fluids. But, there are also some downsides to using low order methods. Increasing the approximation order is more computationally efficient than increasing the mesh density in areas where the solution is expected to be smooth, which for wave simulations is most of the domain away from the free surface. Higher order methods such as discontinuous Galerkin finite element methods (DG FEM) additionally allows more computation to be performed locally, which is beneficial for large distributed calculations [see e.g. @Kubatko09; @Kirby12].

Ocellaris is an exactly incompressible Navier-Stokes solver for free surface flows with a DG FEM based numerical method that supports higher order finite elements and contains specially designed slope limiting stabilisation filters to be able to handle large density transitions [@velslopelim; @OcellarisDoc]. Ocellaris is implemented in Python and C++ with FEniCS [@FEniCSBook; @FEniCS1.5] as the backend for the mesh and finite element assembly. PETSc is used for solving the resulting linear systems [@PETSc].

Ocellaris uses a YAML based input file format documented in the user guide available at [ocellaris.org](https://www.ocellaris.org). The mesh geometry can be defined directly in the input file for simple geometries, or it can be loaded from any file format supported by meshio [@meshio], as long as the unstructured meshes are simplicial; triangles in 2D and tetrahedra in 3D. Due to the flexible nature of the implementation, custom numerical models can be added to a simulation by referencing external Python modules from the input file. Ocellaris uses the XDMF file format [@xdmf] for visualisations and a custom HDF5 format for restart files [@hdf5].

![Cylinder in waves. Rendered in Blender after using Paraview to extract the free surface from an Ocellaris simulation.](https://www.ocellaris.org/figures/cylinder_in_waves.jpg)

# References

