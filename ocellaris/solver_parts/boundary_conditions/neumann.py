# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from . import register_boundary_condition, BoundaryConditionCreator
from ocellaris.utils import CodedExpression, OcellarisCppExpression


DEFAULT_ENFORCE_ZERO_FLUX = False


class OcellarisNeumannBC(object):
    def __init__(
        self, simulation, value, subdomain_id, enforce_zero_flux=DEFAULT_ENFORCE_ZERO_FLUX
    ):
        """
        A simple storage class for Neumann boundary conditions
        """
        self.simulation = simulation
        self._value = value
        self.subdomain_id = subdomain_id
        self.enforce_zero_flux = enforce_zero_flux

    def func(self):
        """
        The boundary value derivative function
        """
        return self._value

    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)

    def __repr__(self):
        return '<OcellarisNeumannBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('ConstantGradient')
class NeumannBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant value Neumann condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Neumann condition
        """
        self.simulation = simulation
        value = inp_dict.get_value('value', required_type='any')
        enforce_zero_flux = inp_dict.get_value(
            'enforce_zero_flux', DEFAULT_ENFORCE_ZERO_FLUX, 'bool'
        )

        if isinstance(value, list):
            assert len(value) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_neumann_condition(name, value[d], subdomain_id)
        else:
            self.register_neumann_condition(var_name, value, subdomain_id, enforce_zero_flux)

    def register_neumann_condition(self, var_name, value, subdomain_id, enforce_zero_flux):
        """
        Add a Neumann condition to this variable
        """
        assert isinstance(value, (float, int))
        df_value = dolfin.Constant(value)

        # Store the boundary condition for use in the solver
        bc = OcellarisNeumannBC(self.simulation, df_value, subdomain_id, enforce_zero_flux)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info('    ConstantGradient %r for %s' % (value, var_name))


@register_boundary_condition('CodedGradient')
class CodedNeumannBoundary(BoundaryConditionCreator):
    description = 'A coded Neumann condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Neumann condition with coded value
        """
        self.simulation = simulation

        # Make a dolfin Expression object that runs the code string
        code = inp_dict.get_value('code', required_type='any')
        enforce_zero_flux = inp_dict.get_value(
            'enforce_zero_flux', DEFAULT_ENFORCE_ZERO_FLUX, 'bool'
        )

        if isinstance(code, list):
            assert len(code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                description = 'coded gradient boundary condition for %s' % name
                sub_code = inp_dict.get_value('code/%d' % d, required_type='string')
                expr = CodedExpression(simulation, sub_code, description)
                self.register_neumann_condition(
                    name, expr, subdomains, subdomain_id, enforce_zero_flux
                )
        else:
            description = 'coded gradient boundary condition for %s' % var_name
            expr = CodedExpression(simulation, code, description)
            self.register_neumann_condition(
                var_name, expr, subdomains, subdomain_id, enforce_zero_flux
            )

    def register_neumann_condition(
        self, var_name, expr, subdomains, subdomain_id, enforce_zero_flux
    ):
        """
        Store the boundary condition for use in the solver
        """
        bc = OcellarisNeumannBC(self.simulation, expr, subdomain_id, enforce_zero_flux)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    Coded gradient for %s' % var_name)


@register_boundary_condition('CppCodedGradient')
class CppCodedNeumannBoundary(BoundaryConditionCreator):
    description = 'A C++ coded Neumann boundary condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        C++ coded Neumann condition
        """
        self.simulation = simulation
        self.func_space = simulation.data['V%s' % var_name]
        cpp_code = inp_dict.get_value('cpp_code', required_type='any')
        enforce_zero_flux = inp_dict.get_value(
            'enforce_zero_flux', DEFAULT_ENFORCE_ZERO_FLUX, 'bool'
        )

        if isinstance(cpp_code, list):
            assert len(cpp_code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                sub_code = inp_dict.get_value('cpp_code/%d' % d, required_type='string')
                self.register_neumann_condition(name, sub_code, subdomain_id, enforce_zero_flux)
        else:
            self.register_neumann_condition(var_name, cpp_code, subdomain_id, enforce_zero_flux)

    def register_neumann_condition(self, var_name, cpp_code, subdomain_id, enforce_zero_flux):
        """
        Add a C++ coded Neumann condition to this variable
        """
        description = 'boundary condititon for %s' % var_name
        P = self.func_space.ufl_element().degree()
        expr = OcellarisCppExpression(self.simulation, cpp_code, description, P, update=True)
        bc = OcellarisNeumannBC(self.simulation, expr, subdomain_id, enforce_zero_flux)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    C++ coded gradient for %s' % var_name)
