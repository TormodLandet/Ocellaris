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

High Reynolds number free surface flows, where the free surface motion does not cause wave breaking, have typically been simulated with potential flow methods. Currently the increasing availability of high performance computing resources and the need to include vorticity and viscous effects is introducing full Navier-Stokes solvers to fields that have traditionally used potential flow methods, such as ship design. Low order numerical methods, such as the finite volume method, are most often used for such Navier-Stokes simulations. Low order finite volume methods can produce exactly incompressible velocity fields, which are important for the advection of the free surface, and does not suffer from Gibbs oscillations near the large jump in density and momentum at the interface between the immiscible fluids.

Ocellaris is an exactly incompressible Navier-Stokes solver for free surface flows with a discontinuous Galerkin finite element based numerical method that supports higher order finite elements and contains specially designed slope limiting stabilisation filters to be able to handle large density transitions [@OcellarisDoc; @OcellarisSrc]. Using high order elements in areas where the solution is expected to be smooth is more computationally efficient than increasing the number of low order computational cells to capture the same smooth field.

Ocellaris is implemented in Python and C++ with FEniCS [@FEniCS] as the backend for the mesh and finite element assembly. PETSc is used for solving the resulting linear systems [@PETSc].

# References
