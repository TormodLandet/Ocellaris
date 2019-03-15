# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from ocellaris.utils import ocellaris_error
from . import register_known_field, KnownField


@register_known_field('SharpField')
class SharpField(KnownField):
    description = 'A field with different values above and below a plane'

    def __init__(self, simulation, field_inp):
        """
        A scalar field
        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a sharp field %r' % self.name)
        simulation.log.info('    Variable: %r' % self.var_name)
        simulation.log.info('    Local proj: %r' % self.local_projection)
        simulation.log.info('    Proj. degree: %r' % self.projection_degree)
        simulation.log.info('    Poly. degree: %r' % self.polynomial_degree)

        mesh = simulation.data['mesh']
        self.V = dolfin.FunctionSpace(mesh, 'DG', self.polynomial_degree)
        self.func = dolfin.Function(self.V)

        if self.local_projection:
            # Represent the jump using a quadrature element
            quad_elem = dolfin.FiniteElement(
                'Quadrature',
                mesh.ufl_cell(),
                self.projection_degree,
                quad_scheme="default",
            )
            xpos, ypos, zpos = self.xpos, self.ypos, self.zpos
            if simulation.ndim == 2:
                cpp0 = 'x[0] < %r and x[1] < %r' % (xpos, ypos)
            else:
                cpp0 = 'x[0] < %r and x[1] < %r and x[2] < %r' % (xpos, ypos, zpos)
            cpp = '(%s) ? %r : %r' % (cpp0, self.val_below, self.val_above)
            e = dolfin.Expression(cpp, element=quad_elem)
            dx = dolfin.dx(metadata={'quadrature_degree': self.projection_degree})

            # Perform the local projection
            u, v = dolfin.TrialFunction(self.V), dolfin.TestFunction(self.V)
            a = u * v * dolfin.dx
            L = e * v * dx
            ls = dolfin.LocalSolver(a, L)
            ls.solve_local_rhs(self.func)

            # Clip values to allowable bounds
            arr = self.func.vector().get_local()
            val_min = min(self.val_below, self.val_above)
            val_max = max(self.val_below, self.val_above)
            clipped_arr = numpy.clip(arr, val_min, val_max)
            self.func.vector().set_local(clipped_arr)
            self.func.vector().apply('insert')

        else:
            # Initialise the sharp static field
            dm = self.V.dofmap()
            arr = self.func.vector().get_local()
            above, below = self.val_above, self.val_below
            xpos, ypos, zpos = self.xpos, self.ypos, self.zpos
            for cell in dolfin.cells(simulation.data['mesh']):
                mp = cell.midpoint()[:]
                dof, = dm.cell_dofs(cell.index())
                if (
                    mp[0] < xpos
                    and mp[1] < ypos
                    and (simulation.ndim == 2 or mp[2] < zpos)
                ):
                    arr[dof] = below
                else:
                    arr[dof] = above
            self.func.vector().set_local(arr)
            self.func.vector().apply('insert')

    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        self.var_name = field_inp.get_value(
            'variable_name', 'phi', required_type='string'
        )
        self.val_above = field_inp.get_value('value_above', required_type='float')
        self.val_below = field_inp.get_value('value_below', required_type='float')
        self.xpos = field_inp.get_value('x', 1e100, required_type='float')
        self.ypos = field_inp.get_value('y', 1e100, required_type='float')
        self.zpos = field_inp.get_value('z', 1e100, required_type='float')
        self.local_projection = field_inp.get_value('local_projection', True, 'bool')
        self.projection_degree = field_inp.get_value('projection_degree', 6, 'bool')
        self.polynomial_degree = field_inp.get_value('polynomial_degree', 0, 'bool')

    def get_variable(self, name):
        if not name == self.var_name:
            ocellaris_error(
                'Sharp field does not define %r' % name,
                'This sharp field defines %r' % self.var_name,
            )
        return self.func
