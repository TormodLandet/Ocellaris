# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from . import register_boundary_condition, BoundaryConditionCreator
from ocellaris.utils import OcellarisError


class OcellarisOutletBC(object):
    def __init__(self, simulation, subdomain_id):
        """
        A simple storage class for zero stress open outlet
        boundary conditions. This is used when defining the
        linear part of the weak forms
        """
        self.simulation = simulation
        self.subdomain_id = subdomain_id
        self.hydrostatic = False
        self.enforce_zero_flux = False

    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)

    def __repr__(self):
        return '<OcellarisOutletBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('OpenOutletBoundary')
class OpenOutletBoundary(BoundaryConditionCreator):
    description = 'An open outlet zero stress boundary condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Open outlet boundary condition
        """
        self.simulation = simulation
        hydrostatic = inp_dict.get_value('hydrostatic', False, 'bool')

        if not var_name == 'u':
            raise OcellarisError(
                'Error in boundary condition',
                'OpenOutletBoundary should be applied to ' '"u" only, p should be left out',
            )

        # Store the boundary condition for use in the solver
        bc = OcellarisOutletBC(self.simulation, subdomain_id)
        bc.hydrostatic = hydrostatic
        self.simulation.data['outlet_bcs'].append(bc)

        self.simulation.log.info('    OpenOutletBoundary for %s' % var_name)
