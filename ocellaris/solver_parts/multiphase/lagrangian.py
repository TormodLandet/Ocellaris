# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from . import register_multi_phase_model, MultiPhaseModel
from .vof import VOFMixin


# Default values, can be changed in the input file
CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION = False


@register_multi_phase_model('Lagrangian')
class LagrangianMeshMorpher(VOFMixin, MultiPhaseModel):
    description = 'A purely Lagrangian multiphase model'

    def __init__(self, simulation):
        """
        A purely Lagrangian multiphase model. The mesh is moved
        according to the calculated fluid velocity after each time
        step. This will obviously distort the mesh in allmost all
        calculations.

        This was implemented as a stepping stone to ALE, and to test
        hydrostatic pressure calculations where the correct answer
        is zero velocity everywhere for all time and ALE should not
        be necessary.

        To initialise the multi phase field the colour function must
        be specified in the input file (as initial condition for "cp").
        The colour function is unity when rho=rho0 and nu=nu0 and
        zero when rho=rho1 and nu=nu1
        """
        self.simulation = simulation

        # Define colour function
        V = simulation.data['Vc']
        c = dolfin.Function(V)
        simulation.data['c'] = c
        simulation.data['cp'] = c  # Initial conditions are specified for cp, so we need this alias

        # Setup the mesh morpher
        simulation.mesh_morpher.setup()

        # Calculate mu from rho and nu (i.e mu is quadratic in c) or directly from c (linear in c)
        self.calculate_mu_directly_from_colour_function = simulation.input.get_value(
            'multiphase_solver/calculate_mu_directly_from_colour_function',
            CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION,
            'bool',
        )

        # Get the physical properties
        self.set_physical_properties(read_input=True)

        # Update the rho and nu fields after each time step
        simulation.hooks.add_post_timestep_hook(self.update, 'LagrangianMeshMorpher - update mesh')

        simulation.log.info('Creating Lagrangian mesh morphing multiphase model')

    def get_colour_function(self, k):
        """
        The colour function follows the cells and does not ever change
        """
        return self.simulation.data['c']

    def update(self):
        """
        Update the mesh position according to the calculated fluid velocities
        """
        timer = dolfin.Timer('Ocellaris Lagrangian mesh update')
        sim = self.simulation

        # Get updated mesh velocity (use the fluid velocity)
        for d in range(sim.ndim):
            ui = sim.data['u%d' % d]
            umi = sim.data['u_mesh%d' % d]
            dolfin.project(ui, V=umi.function_space(), function=umi)

        # Perform the morphing
        sim.mesh_morpher.morph_mesh()

        # Report properties of the colour field
        sum_c = dolfin.assemble(sim.data['c'] * dolfin.dx)
        sim.reporting.report_timestep_value('sum(c)', sum_c)

        timer.stop()
