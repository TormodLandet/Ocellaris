# encoding: utf8
from __future__ import division
import numpy
import dolfin
from solenoidal import SolenoidalLimiter, COST_FUNCTIONS
from ocellaris.utils import verify_key
from ocellaris.solver_parts.boundary_conditions import mark_cell_layers
from ocellaris.solver_parts.slope_limiter import SlopeLimiter
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase


DEFAULT_LIMITER_U = 'HierarchicalTaylor'
DEFAULT_LIMITER_W = 'LocalExtrema'


@register_velocity_slope_limiter('Solenoidal')
class SolenoidalSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'Limit in a divergence free polynomial space'
    
    def __init__(self, simulation, vel, vel_name, vel2, use_cpp=True):
        """
        Use a solenoidal polynomial slope limiter on the convecting velocity field
        w and and use a prelimiter to limit the convected velocity u and set the 
        target for the convecting velocity w
        """
        inp = simulation.input.get_value('slope_limiter/%s' % vel_name, required_type='Input')
        self.additional_plot_funcs = []
        
        # Verify input
        V = vel[0].function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        dim, = vel.ufl_shape
        cost_func = simulation.input.get_value('slope_limiter/cost_function', DEFAULT_LIMITER_W, 'string')
        loc = 'SolenoidalSlopeLimiterVelocity'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (2,), loc)
        verify_key('function dim', dim, (2,), loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        verify_key('cost function', cost_func, COST_FUNCTIONS, loc)
        
        # We expect the convective velocity as vel2 and the convected as vel
        assert vel2 is not None
        
        # Limit all cells regardless of location?
        self.limit_none = inp.get_value('limit_no_cells', False, 'bool')
        
        # Get the IsoSurface probe used to locate the free surface
        self.probe_name = inp.get_value('surface_probe', None, 'string')
        self.surface_probe = None
        self.limit_selected_cells_only = self.probe_name is not None
        
        # Use prelimiter to set (possibly) extended valid bounds and to avoid excessive
        # limiting of the convecting velocity w. The prelimiter is also used for slope 
        # limiting the convected velocity u. This is not treated with the solenoidal limiter
        # since we cannot guarantee that all Gibbs oscillations are removed
        prelimiter_method = inp.get_value('prelimiter', DEFAULT_LIMITER_U, 'string')
        simulation.log.info('        Using prelimiter %r for %s and for solenoidal targets'
                            % (prelimiter_method, vel_name))
        self.prelimiters = []
        for d in range(dim):
            name = '%s%d' % (vel_name, d)
            lim = SlopeLimiter(simulation, vel_name, vel[d], name, method=prelimiter_method)
            self.prelimiters.append(lim)
        
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
        self.sollim = SolenoidalLimiter(V, cost_function=cost_func, use_cpp=use_cpp, cf_options=cf_options)
        self.limit_cell = self.sollim.limit_cell
        
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

        # Boundary cells where we do not fully trust the prelimiter targets (we trust BCs more)
        self.skip_target_cells = mark_cell_layers(simulation, V, layers=1)
        
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
            self.limit_cell[:] = False
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
            self.limit_cell[:Ncells] = surface_cells[:Ncells]

        # Set the solenoidal invariants used in the solenoidal limiting process.
        self.sollim.set_invariants(self.vel)

        # Run prelimiters if they exist and set the target values in vel2
        ui_tmp = self.simulation.data['ui_tmp']
        for i, lim in enumerate(self.prelimiters):
            # Save ui before limiting
            ui_tmp.assign(self.vel[i])
    
            # Limit ui
            lim.limit_cell[:] = self.limit_cell[:]
            lim.run()

            # Set target for optimization of w
            self.vel2[i].assign(self.vel[i])
            
            # Remove boundary cells from targets to avoid problems with BCS
            arr = self.vel2[i].vector().get_local()
            arr_unlim = ui_tmp.vector().get_local()
            V = self.vel2[i].function_space()
            dm = V.dofmap()          
            for cid in self.skip_target_cells:
                for dof in dm.cell_dofs(cid):
                    arr[dof] = arr_unlim[dof]
            self.vel2[i].vector().set_local(arr)
            self.vel2[i].vector().apply('insert')
        
        # Run the optimizing approximate vector slope limiter with vel2 as
        # target and save results to vel2 and well
        self.sollim.run(self.vel2, lim.active, self.max_cost)
        
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
