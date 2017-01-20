# encoding: utf8
from __future__ import division
import numpy
import dolfin
from solenoidal import SolenoidalLimiter
from ocellaris.utils import verify_key
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase


@register_velocity_slope_limiter('Solenoidal')
class SolenoidalSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'Limit in a divergence free polynomial space'
    
    def __init__(self, simulation, vel, vel_name, use_cpp=True):
        """
        Use a solenoidal polynomial slope limiter on the velocity field 
        """
        inp = simulation.input.get_value('slope_limiter/%s' % vel_name, required_type='Input')
        self.additional_plot_funcs = []
        
        # Verify input
        V = vel[0].function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        loc = 'SolenoidalSlopeLimiterVelocity'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (2,), loc)
        verify_key('function shape', vel.ufl_shape, [(2,)], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        
        # Limit all cells regardless of location?
        self.limit_all = inp.get_value('limit_all', False, 'bool')
        
        # Get the IsoSurface probe used to locate the free surface
        self.probe_name = inp.get_value('surface_probe', None, 'string')
        self.surface_probe = self.active_cells = None
        
        # Store input
        self.simulation = simulation
        self.vel = vel
        self.vel_name = vel_name
        self.degree = degree
        self.mesh = mesh
        self.use_cpp = use_cpp
        
        # Create slope limiter
        self.sollim = SolenoidalLimiter(vel)
        
        if self.probe_name is not None or self.limit_all:
            V0 = dolfin.FunctionSpace(self.mesh, 'DG', 0)
            self.active_cells = dolfin.Function(V0)
            aname = 'SolenoidalSlopeLimiterVelocity_%s' % self.vel_name
            self.active_cells.rename(aname, aname)
            self.additional_plot_funcs.append(self.active_cells)
            
            # Cell dofs
            tdim = self.mesh.topology().dim()
            Ncells = self.mesh.topology().ghost_offset(tdim)
            dm0 = V0.dofmap()
            self.cell_dofs_V0 = numpy.array([int(dm0.cell_dofs(i)) for i in xrange(Ncells)], int)
        
        simulation.hooks.add_pre_simulation_hook(self.setup, 'SolenoidalSlopeLimiterVelocity - setup')
    
    def setup(self):
        if self.probe_name is not None:
            verify_key('surface_probe', self.probe_name, self.simulation.probes,
                       'solenoidal slope limiter for %s' % self.vel_name)
            self.surface_probe = self.simulation.probes[self.probe_name]
            self.simulation.log.info('Marking cells for limiting based on probe "%s" for %s'
                                     % (self.surface_probe.name, self.vel_name))
        
        if self.limit_all:
            self.simulation.log.info('Marking all cells for limiting of %s' % self.vel_name)
            self.sollim.limit_cell[:] = True
            self.active_cells.vector()[:] = 1.0
    
    def run(self):
        """
        Perform slope limiting of the velocity field
        """
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
            Ncells = len(self.cell_dofs_V0)
            self.sollim.limit_cell[:Ncells] = surface_cells[:Ncells]
        
        self.sollim.run()
        
        # Update the plot output
        if self.active_cells:
            Ncells = len(self.cell_dofs_V0)
            arr = self.active_cells.vector().get_local()
            arr[self.cell_dofs_V0] = self.sollim.limit_cell[:Ncells]
            self.active_cells.vector().set_local(arr)
            self.active_cells.vector().apply('insert')
