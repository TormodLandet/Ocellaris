# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import ocellaris_error, OcellarisCppExpression
from . import register_known_field, KnownField, DEFAULT_POLYDEG


@register_known_field('ScalarField')
class ScalarField(KnownField):
    description = 'Scalar field'

    def __init__(self, simulation, field_inp):
        """
        A scalar field
        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a scalar field %r' % self.name)
        simulation.log.info('    Variable: %r' % self.var_name)
        simulation.log.info('    Stationary: %r' % self.stationary)
        simulation.log.info('    Polynomial degree: %r' % self.polydeg)

        self.expr = None
        self.func = None
        self.updater = None
        self.V = dolfin.FunctionSpace(simulation.data['mesh'], 'CG', self.polydeg)
        simulation.hooks.add_pre_timestep_hook(
            self.update, 'Update scalar field %r' % self.name
        )

    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        self.var_name = field_inp.get_value(
            'variable_name', 'phi', required_type='string'
        )
        self.stationary = field_inp.get_value('stationary', False, required_type='bool')
        self.polydeg = field_inp.get_value(
            'polynomial_degree', DEFAULT_POLYDEG, required_type='int'
        )
        self.cpp_code = field_inp.get_value('cpp_code', required_type='string')

    def update(self, timestep_number, t, dt):
        """
        Called by simulation.hooks on the start of each time step
        """
        if self.stationary:
            return
        if self.updater is not None:
            self.updater(timestep_number, t, dt)
            self.func.interpolate(self.expr)

    def _get_expression(self):
        if self.expr is None:
            expr, updater = OcellarisCppExpression(
                self.simulation,
                self.cpp_code,
                'Scalar field %r' % self.name,
                self.polydeg,
                update=False,
                return_updater=True,
            )
            self.expr = expr
            self.updater = updater
        return self.expr

    def get_variable(self, name):
        if not name == self.var_name:
            ocellaris_error(
                'Scalar field does not define %r' % name,
                'This scalar field defines %r' % self.var_name,
            )
        if self.func is None:
            expr = self._get_expression()
            func = dolfin.Function(self.V)
            func.interpolate(expr)
            self.func = func
        return self.func
