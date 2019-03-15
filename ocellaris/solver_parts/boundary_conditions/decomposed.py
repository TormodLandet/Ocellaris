# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
Boundary conditions decomposed in wall normal and parallel
directions. Can be used, e.g., for wall slip on walls that
are not parallel to a cartesian axis
"""
from ocellaris.utils import OcellarisError
from . import register_boundary_condition, BoundaryConditionCreator


class FreeSlipBC(object):
    def __init__(self, simulation, subdomain_id):
        """
        A simple storage class for free slip boundaray conditions.
        """
        self.simulation = simulation
        self.subdomain_id = subdomain_id

    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)

    def __repr__(self):
        return '<OcellarisFreeSlipBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('FreeSlip')
class FreeSlipBoundary(BoundaryConditionCreator):
    description = 'Free slip on non-penetrable wall'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Free slip boundary condition
        """
        if var_name[-1] in '0123456789':
            raise OcellarisError(
                'Error in FreeSlip boundary conditions',
                'You must give the name of a vector field, like "u", '
                'not a component. You gave %r' % var_name,
            )

        # Store the boundary condition for use in the solver
        bc = FreeSlipBC(simulation, subdomain_id)
        bcs = simulation.data['slip_bcs']
        bcs.setdefault(var_name, []).append(bc)

        simulation.log.info('    FreeSlipBC for %s' % var_name)
