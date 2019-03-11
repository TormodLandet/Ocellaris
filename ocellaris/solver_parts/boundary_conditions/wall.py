# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from . import register_boundary_condition, BoundaryConditionCreator
from .neumann import OcellarisNeumannBC


@register_boundary_condition('WallPressure')
class WallPressureBoundaryCondition(BoundaryConditionCreator):
    description = 'Boundary condition for pressure at a wall, dp/dn = ρ g⋅n'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Wall pressure boundary condition
        """
        mesh = simulation.data['mesh']
        n = dolfin.FacetNormal(mesh)
        g = simulation.data['g']
        rho = simulation.data['rho']
        bc_val = rho * dolfin.inner(n, g)

        # Store the boundary condition for use in the solver
        bc = OcellarisNeumannBC(simulation, bc_val, subdomain_id)
        bcs = simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)

        simulation.log.info('    WallPressure for %s' % var_name)
