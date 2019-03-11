# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import ocellaris_error, OcellarisCppExpression
from . import register_known_field, KnownField
from .base_wave_field import COLOUR_PROJECTION_DEGREE


@register_known_field('WaterBlock')
class WaterBlockField(KnownField):
    description = 'A block of water'

    def __init__(self, simulation, field_inp):
        """
        A block of water, projected to obtain the best representation of
        the colour field c given a mesh that does not conform to the block
        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a water block field %r' % self.name)
        simulation.log.info('    Variable: %r' % self.var_name)
        simulation.log.info('    Polynomial degree: %r' % self.polydeg)
        simulation.log.info('    C++ code: %s' % self.cpp_code)
        simulation.log.info('    Quadrature degree: %r' % self.colour_projection_degree)

        self.func = None
        fs = 'DG' if self.polydeg == 0 else 'CG'
        self.V = dolfin.FunctionSpace(simulation.data['mesh'], fs, self.polydeg)

    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        self.var_name = field_inp.get_value('variable_name', 'c', 'string')
        self.polydeg = field_inp.get_value('polynomial_degree', 0, 'int')

        self.xmin = field_inp.get_value('xmin', 0, 'float')
        self.xmax = field_inp.get_value('xmax', 0, 'float')
        self.ymin = field_inp.get_value('ymin', 0, 'float')
        self.ymax = field_inp.get_value('ymax', 0, 'float')
        self.zmin = field_inp.get_value('zmin', 0, 'float')
        self.zmax = field_inp.get_value('zmax', 0, 'float')

        # Project the colour function to DG0 (set degree to -1 to prevent this)
        self.colour_projection_degree = field_inp.get_value(
            'colour_projection_degree', COLOUR_PROJECTION_DEGREE, 'int'
        )

        # Construct C++ code
        cpp = []
        if self.xmin != self.xmax:
            cpp.append('x[0] >= %r' % self.xmin)
            cpp.append('x[0] <= %r' % self.xmax)
        if self.ymin != self.ymax:
            cpp.append('x[1] >= %r' % self.ymin)
            cpp.append('x[1] <= %r' % self.ymax)
        if self.zmin != self.zmax:
            cpp.append('x[2] >= %r' % self.zmin)
            cpp.append('x[2] <= %r' % self.zmax)
        cpp = ['(%s)' % part for part in cpp]
        self.cpp_code = '(%s) ? 1.0 : 0.0' % ' and '.join(cpp)

    def get_variable(self, name):
        if not name == self.var_name:
            ocellaris_error(
                'Scalar field does not define %r' % name,
                'This scalar field defines %r' % self.var_name,
            )

        if self.func is not None:
            return self.func

        if self.colour_projection_degree >= 0:
            quad_degree = self.colour_projection_degree
            degree = None
        else:
            quad_degree = None
            degree = self.polydeg

        expr = OcellarisCppExpression(
            self.simulation,
            self.cpp_code,
            'Scalar field %r' % self.name,
            degree=degree,
            update=False,
            quad_degree=quad_degree,
        )
        self.func = dolfin.Function(self.V)

        if quad_degree is not None:

            # Ensure that we can use the DG0 trick of dividing by the mass
            # matrix diagonal to get the projected value
            if self.V.ufl_element().degree() != 0:
                ocellaris_error(
                    'Error in wave field %r input' % self.name,
                    'The colour_to_dg_degree projection is '
                    'currently only implemented when c is DG0',
                )
            v = dolfin.TestFunction(self.V)
            dx = dolfin.dx(metadata={'quadrature_degree': quad_degree})
            d = dolfin.CellVolume(self.V.mesh())  # mass matrix diagonal
            form = expr * v / d * dx
            compiled_form = dolfin.Form(form)

            # Perform projection by assembly (DG0 only!)
            dolfin.assemble(compiled_form, tensor=self.func.vector())
        else:
            # Perform standard interpolation
            self.func.interpolate(expr)

        return self.func
