# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from ocellaris.utils import OcellarisError, RunnablePythonString
from . import register_multi_phase_model, MultiPhaseModel
from .vof import VOFMixin


@register_multi_phase_model('HeightFunction')
class HeightFunctionMultiphaseModel(VOFMixin, MultiPhaseModel):
    description = 'Separate the phases with a height function inside a static mesh'
    calculate_mu_directly_from_colour_function = False
    default_polynomial_degree_colour = 1

    def __init__(self, simulation):
        """
        A simple height function multiphase model
        """
        self.simulation = simulation
        simulation.log.info('Creating height function multiphase model')
        assert simulation.ncpu == 1, 'HeightFunctionMultiphaseModel does not run in parallel'

        # Define colour function
        V = simulation.data['Vc']
        simulation.data['c'] = dolfin.Function(V)

        # The free surface mesh
        xmin = simulation.input.get_value(
            'multiphase_solver/height_function_xmin', required_type='float'
        )
        xmax = simulation.input.get_value(
            'multiphase_solver/height_function_xmax', required_type='float'
        )
        Nx = simulation.input.get_value(
            'multiphase_solver/height_function_Nx', required_type='float'
        )
        hinit = simulation.input.get_value(
            'multiphase_solver/height_initial_h', required_type='float'
        )
        self.eps = simulation.input.get_value(
            'multiphase_solver/surface_thickness', required_type='float'
        )
        hmin_code = simulation.input.get_value(
            'multiphase_solver/height_min_code', required_type='string'
        )
        hmax_code = simulation.input.get_value(
            'multiphase_solver/height_max_code', required_type='string'
        )

        # Runnable code, must define 'hmin' or 'hmax' arrays given input xpos array
        self.hmin = RunnablePythonString(
            simulation, hmin_code, 'multiphase_solver/height_min_code', 'hmin'
        )
        self.hmax = RunnablePythonString(
            simulation, hmax_code, 'multiphase_solver/height_max_code', 'hmax'
        )

        # Store colour function dof coordinates
        mesh = simulation.data['mesh']
        gdim = mesh.geometry().dim()
        self.dof_pos = V.tabulate_dof_coordinates().reshape((-1, gdim))

        # Create dummy dolfin Function to hold the height function (for restart file support)
        if Nx > V.dim():
            raise OcellarisError(
                'height_function_Nx too large',
                'Maximum Nx for this colour function is %d' % V.dim(),
            )
        simulation.data['height_function'] = hfunc = dolfin.Function(V)
        hfunc.vector()[:] = hinit
        hfunc.vector().apply('insert')
        simulation.data['height_function_x'] = numpy.linspace(xmin, xmax, Nx)
        self._update_colour()
        simulation.log.info('    Created surface mesh with %d elements' % (Nx - 1))

        # Get the physical properties
        self.set_physical_properties(read_input=True)

        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(
            self.update, 'HeightFunctionMultiphaseModel.update()'
        )
        simulation.hooks.register_custom_hook_point('MultiPhaseModelUpdated')

    def get_colour_function(self, k):
        """
        The colour function follows the cells and does not ever change
        """
        return self.simulation.data['c']

    def _update_colour(self):
        """
        Update the colour function based on the current heigh function
        """
        sim = self.simulation
        hfunc = sim.data['height_function']
        xpos = sim.data['height_function_x']
        Nx = xpos.size
        ypos = hfunc.vector().get_local()[:Nx]

        # Find height at each dof x position
        dofs_x = self.dof_pos[:, 1]
        dofs_y = self.dof_pos[:, 1]
        fs_dist = dofs_y - numpy.interp(dofs_x, xpos, ypos)

        # Update colour function based on dof y position
        c = numpy.tanh(-fs_dist * 4 / self.eps) * 0.5 + 0.5
        sim.data['c'].vector().set_local(c)
        sim.data['c'].vector().apply('insert')

    def update(self, timestep_number, t, dt):
        """
        Update the mesh position according to the calculated fluid velocities
        """
        sim = self.simulation
        hfunc = sim.data['height_function']
        xpos = sim.data['height_function_x']
        Nx = xpos.size
        old_h = hfunc.vector().get_local()

        with dolfin.Timer('Ocellaris update height and colour'):
            # Get updated mesh velocity in each xpos
            all_vel = numpy.zeros_like(xpos)
            vel = numpy.zeros(1, float)
            pos = numpy.zeros(2, float)
            uvert = sim.data['u1']
            for i, x in enumerate(xpos):
                pos[:] = (x, old_h[i])
                try:
                    uvert.eval(vel, pos)
                except Exception:
                    raise OcellarisError(
                        'Error in u1 eval in height function',
                        'Position %r outside the domain?' % (pos,),
                    )
                all_vel[i] = vel[0]

            # Get height function max and min
            hmin = self.hmin.run(xpos=xpos)
            hmax = self.hmax.run(xpos=xpos)
            assert hmin.shape == xpos.shape
            assert hmax.shape == xpos.shape

            # Update the height function
            dt = sim.input.get_value('time/dt', required_type='float')
            old_h[:Nx] += dt * all_vel
            old_h[:Nx] = numpy.clip(old_h[:Nx], hmin, hmax)

            hfunc.vector().set_local(old_h)
            hfunc.vector().apply('insert')

            self._update_colour()

        # Report on the new height function
        self.simulation.reporting.report_timestep_value('min(h)', old_h[:Nx].min())
        self.simulation.reporting.report_timestep_value('max(h)', old_h[:Nx].max())
        self.simulation.hooks.run_custom_hook('MultiPhaseModelUpdated')
