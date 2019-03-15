# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from . import register_boundary_condition, BoundaryConditionCreator


DEFAULT_BLEND = 1.0
DEFAULT_DVAL = 0.0
DEFAULT_NVAL = 0.0
DEFAULT_NSCALE = 1.0
DEFAULT_ENFORCE_ZERO_FLUX = False


class OcellarisRobinBC:
    def __init__(
        self,
        simulation,
        V,
        blend,
        dirichlet_value,
        neumann_value,
        subdomain_marker,
        subdomain_id,
        updater=None,
    ):
        """
        A storage class for a Robin boundary conditions on the form

          n⋅∇φ = 1/b (φ0 - φ) + g

        Where b is a blending parameter and u0 and g are the Dirichlet
        and Neumann values for b → 0 and b → ∞ respectively.
        """
        self.simulation = simulation
        self._blend = blend
        self._dfunc = dirichlet_value
        self._nfunc = neumann_value
        self.subdomain_marker = subdomain_marker
        self.subdomain_id = subdomain_id
        self.enforce_zero_flux = DEFAULT_ENFORCE_ZERO_FLUX
        self._updater = updater

    def blend(self):
        return self._blend

    def dfunc(self):
        return self._dfunc

    def nfunc(self):
        return self._nfunc

    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)

    def update(self):
        """
        Update the time and other parameters used in the BC.
        This is used every timestep and for all RK substeps
        """
        if self._updater:
            self._updater(self.simulation.timestep, self.simulation.time, self.simulation.dt)

    def __repr__(self):
        return '<OcellarisRobinBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('ConstantRobin')
class ConstantRobinBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant value Robin boundary condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Robin condition with constant values
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]

        blend = inp_dict.get_value('blend', DEFAULT_BLEND, required_type='any')
        dval = inp_dict.get_value('dval', DEFAULT_DVAL, required_type='any')
        nval = inp_dict.get_value('nval', DEFAULT_NVAL, required_type='any')
        if isinstance(blend, list):
            assert len(blend) == simulation.ndim
            assert len(dval) == simulation.ndim
            assert len(nval) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_robin_condition(
                    name, blend[d], dval[d], nval[d], subdomains, subdomain_id
                )
        else:
            self.register_robin_condition(var_name, blend, dval, nval, subdomains, subdomain_id)

    def register_robin_condition(self, var_name, blend, dval, nval, subdomains, subdomain_id):
        """
        Add a Dirichlet condition to this variable
        """
        assert isinstance(blend, (float, int))
        assert isinstance(dval, (float, int))
        assert isinstance(nval, (float, int))
        df_blend = dolfin.Constant(blend)
        df_dval = dolfin.Constant(dval)
        df_nval = dolfin.Constant(nval)

        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(
            self.simulation, self.func_space, df_blend, df_dval, df_nval, subdomains, subdomain_id
        )
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info(
            '    Constant Robin BC "dX/dn = 1/%r (%r - X) + %r" for X = %s'
            % (blend, dval, nval, var_name)
        )
