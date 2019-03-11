# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import OcellarisError
from . import register_known_field, KnownField


@register_known_field('FreeSurfaceZone')
class FreeSurfaceZone(KnownField):
    description = '1.0 near the interface and 0.0 away from the interface'

    def __init__(self, simulation, field_inp):
        """
        Mark the area within a distance radius from the free surface
        """
        self.simulation = simulation
        self._construct(field_inp)

    def _construct(self, field_inp):
        # Read the input
        self.name = field_inp.get_value('name', required_type='string')
        self.var_name = field_inp.get_value('variable_name', 'phi', 'string')
        self.stationary = field_inp.get_value('stationary', False, 'bool')
        self.radius = field_inp.get_value('radius', required_type='float')
        self.plot = field_inp.get_value('plot', False, 'bool')

        # Show the input
        sim = self.simulation
        sim.log.info('Creating a free surface zone field %r' % self.name)
        sim.log.info('    Variable: %r' % self.var_name)
        sim.log.info('    Stationary: %r' % self.stationary)

        # Create the level set model
        mesh = sim.data['mesh']
        func_name = '%s_%s' % (self.name, self.var_name)
        self.V = dolfin.FunctionSpace(mesh, 'DG', 0)
        self.function = dolfin.Function(self.V)
        self.function.rename(func_name, func_name)
        sim.data[func_name] = self.function

        # Get the level set view
        level_set_view = sim.multi_phase_model.get_level_set_view()
        level_set_view.add_update_callback(self.update)

        # Form to compute the cell average distance to the free surface
        v = dolfin.TestFunction(self.V)
        ls = level_set_view.level_set_function
        cv = dolfin.CellVolume(self.V.mesh())
        self.dist_form = dolfin.Form(ls * v / cv * dolfin.dx)

        if self.plot:
            sim.io.add_extra_output_function(self.function)

    def get_variable(self, name):
        if not name == self.var_name:
            raise OcellarisError(
                'FreeSurfaceZone field does not define %r' % name,
                'This FreeSurfaceZone field defines %r' % self.var_name,
            )
        return self.function

    def update(self, force=False):
        # No need to update stationary fields after first update
        if self.stationary and not force:
            return

        # Compute the cell average distance to the free surface
        f = self.function
        dolfin.assemble(self.dist_form, tensor=f.vector())

        # Find the inside cells and the cells in the smoothing zone
        r = abs(f.vector().get_local()) / self.radius
        inside = r <= 1.0
        outside = r > 2.0

        # Create the inside and smoothing zone with a smooth transition
        # from 1 to 0 when r goes from 1 to 2. The slope in both ends is 0
        phi = 2 * r ** 3 - 9 * r ** 2 + 12 * r - 4
        phi[inside] = 1.0
        phi[outside] = 0.0

        f.vector().set_local(phi)
        f.vector().apply('insert')
