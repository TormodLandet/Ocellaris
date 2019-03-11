# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
from ocellaris.utils import ocellaris_error, verify_key
from ocellaris.solver_parts.slope_limiter import SlopeLimiter
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase
from .velocity_limiter_helpers import create_component_limiters


@register_velocity_slope_limiter('Componentwise')
class ComponentwiseSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'Scalar limiting of each component'

    def __init__(self, simulation, vel_u, vel_name, vel_w=None, use_cpp=True):
        """
        Use a standard slope limiter on each component of the velocity field
        """
        self.simulation = simulation
        self.vel_name = vel_name
        inp = simulation.input.get_value('slope_limiter/%s' % vel_name, {}, 'Input')
        comp_method = inp.get_value('comp_method', required_type='string')

        self.vel_u = vel_u
        self.vel_w = vel_w
        self.limit_conv = inp.get_value('limit_conv', False, 'bool')

        # Get the IsoSurface probe used to locate the free surface
        self.probe_name = inp.get_value('surface_probe', None, 'string')
        self.surface_probe = None
        self.limit_selected_cells_only = self.probe_name is not None

        # Create slope limiters
        dim, = vel_u.ufl_shape
        self.limiters = create_component_limiters(
            simulation, vel_name, vel_u, comp_method
        )

        # Check that we can limit only certain cells
        if self.limit_selected_cells_only:
            for lim in self.limiters:
                if not hasattr(lim, 'limit_cell'):
                    ocellaris_error(
                        'Cannot limit only selected cells for %s' % vel_name,
                        'Limiter %r does not support limiting only selected cells'
                        % type(lim).__name__,
                    )

        simulation.hooks.add_pre_simulation_hook(
            self.setup, 'ComponentwiseSlopeLimiterVelocity - setup'
        )

    def setup(self):
        """
        Deferred setup tasks that are run after the Navier-Stokes solver has finished its setup
        """
        if self.probe_name is not None:
            verify_key(
                'surface_probe',
                self.probe_name,
                self.simulation.probes,
                'componentwise slope limiter for %s' % self.vel_name,
            )
            self.surface_probe = self.simulation.probes[self.probe_name]
            self.simulation.log.info(
                'Marking cells for limiting based on probe "%s" for %s'
                % (self.surface_probe.name, self.vel_name)
            )

    def run(self):
        """
        Perform slope limiting of the velocity field
        """
        # Check if we should limit only cells near an IsoSurface
        if self.surface_probe is not None:
            connectivity_CF = self.simulation.data['connectivity_CF']
            connectivity_FC = self.simulation.data['connectivity_FC']
            surface_cells = self.surface_probe.cells_with_surface.copy()
            active_cells = numpy.zeros_like(surface_cells)

            # Mark neighbours of the surface cells in Nlayers layers
            Nlayers = 2
            for _ in range(Nlayers):
                for cid, active in enumerate(surface_cells):
                    if not active:
                        continue

                    for fid in connectivity_CF(cid):
                        for nid in connectivity_FC(fid):
                            active_cells[nid] = True
                surface_cells[:] = active_cells

            # Mark cells to be limited
            for lim in self.limiters:
                Ncells = len(lim.limit_cell)
                lim.limit_cell[:] = surface_cells[:Ncells]

        # Perform limiting
        for lim in self.limiters:
            lim.run()

        # Apply limiting also to the convective field?
        if self.limit_conv:
            dim, = self.vel_u.ufl_shape
            for d in range(dim):
                self.vel_w[d].assign(self.vel_u[d])
