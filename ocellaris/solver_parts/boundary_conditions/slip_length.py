# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import OcellarisCppExpression, OcellarisError, verify_field_variable_definition
from .robin import OcellarisRobinBC
from . import register_boundary_condition, BoundaryConditionCreator


def df_wrap(val, description, degree, sim):
    """
    Wrap numbers as dolfin.Constant and strings as
    dolfin.Expression C++ code. Lists must be ndim
    long and contain either numbers or strings
    """
    if isinstance(val, (int, float)):
        # A real number
        return dolfin.Constant(val)
    elif isinstance(val, str):
        # A C++ code string
        return OcellarisCppExpression(sim, val, description, degree)
    elif isinstance(val, (list, tuple)):
        D = sim.ndim
        L = len(val)
        if L != D:
            raise OcellarisError(
                'Invalid length of list',
                'BC list in "%r" must be length %d, is %d.' % (description, D, L),
            )

        if all(isinstance(v, str) for v in val):
            # A list of C++ code strings
            return OcellarisCppExpression(sim, val, description, degree)
        else:
            # A mix of constants and (possibly) C++ strings
            val = [df_wrap(v, description + ' item %d' % i, degree, sim) for i, v in enumerate(val)]
            return dolfin.as_vector(val)


@register_boundary_condition('SlipLength')
class SlipLengthBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant slip length (Navier) boundary condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Wall slip length (Navier) boundary condition with constant value
        """
        self.simulation = simulation
        vn = var_name[:-1] if var_name[-1].isdigit() else var_name
        self.func_space = simulation.data['V%s' % vn]
        dim = self.func_space.num_sub_spaces()
        default_base = 0.0 if dim == 0 else [0.0] * dim

        length = inp_dict.get_value('slip_length', required_type='any')
        base = inp_dict.get_value('value', default_base, required_type='any')
        self.register_slip_length_condition(var_name, length, base, subdomains, subdomain_id)

    def register_slip_length_condition(self, var_name, length, base, subdomains, subdomain_id):
        """
        Add a Robin boundary condition to this variable
        """
        degree = self.func_space.ufl_element().degree()
        df_blend = df_wrap(length, 'slip length for %s' % var_name, degree, self.simulation)
        df_dval = df_wrap(base, 'boundary condition for %s' % var_name, degree, self.simulation)
        df_nval = 0.0

        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(
            self.simulation, self.func_space, df_blend, df_dval, df_nval, subdomains, subdomain_id
        )
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info(
            '    Constant slip length = %r (base %r) for %s' % (length, base, var_name)
        )


@register_boundary_condition('InterfaceSlipLength')
class InterfaceSlipLengthBoundary(BoundaryConditionCreator):
    description = 'A variable slip length (Navier) boundary condition changing along the boundary'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Wall slip length (Navier) boundary condition where the slip length is multiplied
        by a slip factor ∈ [0, 1] that varies along the domain boundary depending on the
        distance to an interface (typically a free surface between two fluids).
        """
        self.simulation = simulation
        vn = var_name[:-1] if var_name[-1].isdigit() else var_name
        self.func_space = simulation.data['V%s' % vn]
        dim = self.func_space.num_sub_spaces()
        default_base = 0.0 if dim == 0 else [0.0] * dim

        length = inp_dict.get_value('slip_length', required_type='any')
        base = inp_dict.get_value('value', default_base, 'float')
        vardef = inp_dict.get_value('slip_factor_function', required_type='string')

        # Use a subdomain as the slip factor (1.0 inside the domain and
        # 0.0 outside with a smooth transition)
        fac = verify_field_variable_definition(simulation, vardef, 'InterfaceSlipLength')
        fac_name = 'subdomain %s' % vardef.split('/')[0]

        self.register_slip_length_condition(
            var_name, length, fac, fac_name, base, subdomains, subdomain_id
        )

    def register_slip_length_condition(
        self, var_name, length, factor, factor_name, base, subdomains, subdomain_id
    ):
        """
        Add a Robin boundary condition to this variable
        """
        degree = self.func_space.ufl_element().degree()
        df_length = df_wrap(length, 'slip length for %s' % var_name, degree, self.simulation)
        df_blend = df_length * factor
        df_dval = df_wrap(base, 'boundary condition for %s' % var_name, degree, self.simulation)
        df_nval = 0.0

        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(
            self.simulation, self.func_space, df_blend, df_dval, df_nval, subdomains, subdomain_id
        )
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info(
            '    Variable slip length %r (base %r) for %s' % (factor_name, base, var_name)
        )

