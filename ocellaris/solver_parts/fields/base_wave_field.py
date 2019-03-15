# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from collections import OrderedDict
import dolfin
from ocellaris.utils import ocellaris_error, OcellarisCppExpression, verify_key
from . import KnownField, DEFAULT_POLYDEG


# Degree used to go from colour field with jump to DG colour field
# Set it to -1 to use the same interpolation as for the other fields
COLOUR_PROJECTION_DEGREE = 6


class BaseWaveField(KnownField):
    def __init__(self, simulation, field_inp):
        """
        This field is not used directly, but as a basis for the Airy and Raschii
        wave fields since they have a lot of functionality in common
        """
        self.simulation = simulation
        self._dependent_fields = []
        self._cpp = {}
        self._expressions = OrderedDict()
        self._functions = OrderedDict()

        self.name = field_inp.get_value('name', required_type='string')
        self.polydeg = field_inp.get_value(
            'polynomial_degree', DEFAULT_POLYDEG, required_type='int'
        )
        self.V = dolfin.FunctionSpace(simulation.data['mesh'], 'CG', self.polydeg)
        simulation.hooks.add_pre_timestep_hook(
            self.update, 'Update wave field %r' % self.name
        )

    def update(self, timestep_number, t, dt):
        """
        Called by simulation.hooks on the start of each time step
        """
        if self.stationary:
            return

        # Update C++ expressions
        for name, func in self._functions.items():
            expr, updater = self._expressions[name]
            updater(timestep_number, t, dt)
            self._interp(name, expr, func)

        # Update dependent fields
        for f in self._dependent_fields:
            f.update(timestep_number, t, dt)

    def register_dependent_field(self, field):
        """
        We may have dependent fields, typically still water outflow
        fields that use our functions in order to compute their own
        functions
        """
        self._dependent_fields.append(field)

    def _get_expression(self, name, quad_degree=None):
        keys = list(self._cpp.keys()) + ['u']
        verify_key('variable', name, keys, 'Raschii wave field %r' % self.name)
        if name not in self._expressions:
            degree = self.polydeg if quad_degree is None else None
            expr, updater = OcellarisCppExpression(
                self.simulation,
                self._cpp[name],
                'Raschii wave field %r' % name,
                degree,
                update=False,
                return_updater=True,
                quad_degree=quad_degree,
            )
            self._expressions[name] = expr, updater
        return self._expressions[name][0]

    def _interp(self, name, expr=None, func=None):
        """
        Interpolate the expression into the given function
        """
        sim = self.simulation

        # Determine the function space
        V = self.V
        quad_degree = None
        if name == 'c' and self.colour_projection_degree >= 0:
            # Perform projection for c instead of interpolation to better
            # capture the discontinuous nature of the colour field
            if 'Vc' not in sim.data:
                ocellaris_error(
                    'Error in wave field %r input' % self.name,
                    'Cannot specify colour_to_dg_degree when c is '
                    'not used in the multiphase_solver.',
                )
            V = sim.data['Vc']
            quad_degree = self.colour_projection_degree

        # Get the expression
        if expr is None:
            expr = self._get_expression(name, quad_degree)

        # Get the function
        if func is None:
            if name not in self._functions:
                self._functions[name] = dolfin.Function(V)
            func = self._functions[name]

        if quad_degree is not None:
            if self.colour_projection_form is None:
                # Ensure that we can use the DG0 trick of dividing by the mass
                # matrix diagonal to get the projected value
                if V.ufl_element().degree() != 0:
                    ocellaris_error(
                        'Error in wave field %r input' % self.name,
                        'The colour_to_dg_degree projection is '
                        'currently only implemented when c is DG0',
                    )
                v = dolfin.TestFunction(V)
                dx = dolfin.dx(metadata={'quadrature_degree': quad_degree})
                d = dolfin.CellVolume(V.mesh())  # mass matrix diagonal
                form = expr * v / d * dx
                self.colour_projection_form = dolfin.Form(form)

            # Perform projection by assembly (DG0 only!)
            dolfin.assemble(self.colour_projection_form, tensor=func.vector())
        else:
            # Perform standard interpolation
            func.interpolate(expr)

    def get_variable(self, name):
        """
        Return a dolfin Function or as_vector[Function ...]) representing
        the given name
        """
        if name == 'u':
            # Assume uhoriz == u0 and uvert = uD where D is 2 or 3
            if self.simulation.ndim == 2:
                return dolfin.as_vector(
                    [self.get_variable('uhoriz'), self.get_variable('uvert')]
                )
            else:
                return dolfin.as_vector(
                    [
                        self.get_variable('uhoriz'),
                        dolfin.Constant(0.0),
                        self.get_variable('uvert'),
                    ]
                )

        if name not in self._functions:
            self._interp(name)
        return self._functions[name]
