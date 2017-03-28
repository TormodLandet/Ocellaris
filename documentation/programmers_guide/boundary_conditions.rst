Boundary conditions
===================

The boundary condition code will both identify regions of the boundary given by
the user in the input file and create boundary condition objects for each 
function (velocity, pressure ...) in this region.

.. autoclass:: ocellaris.solver_parts.boundary_conditions.BoundaryRegion

.. autoclass:: ocellaris.solver_parts.boundary_conditions.BoundaryCondition
