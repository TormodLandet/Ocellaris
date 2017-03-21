# encoding: utf8
from __future__ import division
import numpy
import dolfin
from solenoidal import SolenoidalLimiter, COST_FUNCTIONS
from ocellaris.utils import verify_key
from ocellaris.solver_parts.slope_limiter import SlopeLimiter
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase


@register_velocity_slope_limiter('Solenoidal')
class SolenoidalSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'Limit in a divergence free polynomial space'
    
    def __init__(self, simulation, vel, vel_name, vel2=None, use_cpp=True):
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
        cost_func = simulation.input.get_value('slope_limiter/cost_function', 'LocalExtrema', 'string')
        loc = 'SolenoidalSlopeLimiterVelocity'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (2,), loc)
        verify_key('function shape', vel.ufl_shape, [(2,)], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        verify_key('cost function', cost_func, COST_FUNCTIONS, loc)
        
        # Limit all cells regardless of location?
        self.limit_none = inp.get_value('limit_no_cells', False, 'bool')
        
        # Get the IsoSurface probe used to locate the free surface
        self.probe_name = inp.get_value('surface_probe', None, 'string')
        self.surface_probe = None
        self.limit_selected_cells_only = self.probe_name is not None
        
        # Use prelimiter to set (possibly) extended valid bound and avoid excessive limiting
        prelimiter_method = inp.get_value('prelimiter', None, 'string')
        prelimiter = None
        if prelimiter_method:
            def run_prelim():
                for lim in prelimiters:
                    lim.limit_cell[:] = self.sollim.limit_cell[:]
                    lim.run()
            prelimiters = []
            dim, = vel.ufl_shape
            for d in range(dim):
                name = '%s%d' % (vel_name, d)
                lim = SlopeLimiter(simulation, vel_name, vel[d], name, method=prelimiter_method)
                prelimiters.append(lim)
            prelimiter = run_prelim
        
        # Cost function options
        cf_options = {}
        cf_option_keys = ('max_cost', 'out_of_bounds_penalty_fac', 'out_of_bounds_penalty_const')
        for key in cf_option_keys:
            if key in inp:
                cf_options[key] = inp.get_value(key, required_type='float')
        max_cost = cf_options.pop('max_cost', None)
        
        # Store input
        self.simulation = simulation
        self.vel = vel
        self.vel2 = vel2
        self.vel_name = vel_name
        self.degree = degree
        self.mesh = mesh
        self.use_cpp = use_cpp
        self.cf_options = cf_options
        self.max_cost = max_cost
        
        # Create slope limiter
        self.sollim = SolenoidalLimiter(vel, cost_function=cost_func, use_cpp=use_cpp,
                                        prelimiter=prelimiter, cf_options=cf_options,
                                        max_cost=max_cost, vec_field2=vel2)
        
        # Create plot output functions
        V0 = dolfin.FunctionSpace(self.mesh, 'DG', 0)
        
        # Final cost from minimizer per cell
        self.cell_cost = dolfin.Function(V0)
        cname = 'SolenoidalCostFunc_%s' % self.vel_name
        self.cell_cost.rename(cname, cname)
        self.additional_plot_funcs.append(self.cell_cost)
        
        self.active_cells = None
        if self.limit_selected_cells_only:
            self.active_cells = dolfin.Function(V0)
            aname = 'SolenoidalActiveCells_%s' % self.vel_name
            self.active_cells.rename(aname, aname)
            self.additional_plot_funcs.append(self.active_cells)
        
        # Cell dofs
        tdim = self.mesh.topology().dim()
        Ncells = self.mesh.topology().ghost_offset(tdim)
        dm0 = V0.dofmap()
        self.cell_dofs_V0 = numpy.array([int(dm0.cell_dofs(i)) for i in xrange(Ncells)], int)
        
        simulation.hooks.add_pre_simulation_hook(self.setup, 'SolenoidalSlopeLimiterVelocity - setup')
    
    def setup(self):
        """
        Deferred setup tasks that are run after the Navier-Stokes solver has finished its setup
        """
        if self.probe_name is not None:
            verify_key('surface_probe', self.probe_name, self.simulation.probes,
                       'solenoidal slope limiter for %s' % self.vel_name)
            self.surface_probe = self.simulation.probes[self.probe_name]
            self.simulation.log.info('Marking cells for limiting based on probe "%s" for %s'
                                     % (self.surface_probe.name, self.vel_name))
        
        if self.limit_none:
            self.simulation.log.info('Marking no cells for limiting of %s' % self.vel_name)
            self.sollim.limit_cell[:] = False
            self.surface_probe = None
            self.limit_selected_cells_only = False
    
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
        
        # Update plot of cost function
        Ncells = len(self.cell_dofs_V0)
        arr = self.cell_cost.vector().get_local()
        arr[self.cell_dofs_V0] = self.sollim.cell_cost[:Ncells]
        self.cell_cost.vector().set_local(arr)
        self.cell_cost.vector().apply('insert')
        
        # Update the plot output of active cells
        if self.limit_selected_cells_only:
            arr = self.active_cells.vector().get_local()
            arr[self.cell_dofs_V0] = self.sollim.limit_cell[:Ncells]
            self.active_cells.vector().set_local(arr)
            self.active_cells.vector().apply('insert')
