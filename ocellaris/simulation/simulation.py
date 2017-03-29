import time
import numpy
import dolfin
from ocellaris.utils import ocellaris_error, velocity_change
from ocellaris.utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data
from .hooks import Hooks
from .input import Input
from .plotting import Plotting
from .reporting import Reporting
from .log import Log
from .io import InputOutputHandling
from .solution_properties import SolutionProperties
from .setup import setup_simulation


class Simulation(object):
    def __init__(self):
        """
        Represents one Ocellaris simulation. The Simulation class 
        connects the input file, geometry, mesh and more with the
        solver and the result plotting and reporting tools     
        """
        self.ncpu = dolfin.MPI.size(dolfin.mpi_comm_world())
        self.rank = dolfin.MPI.rank(dolfin.mpi_comm_world())
        
        self.hooks = Hooks(self)
        self.input = Input(self)
        self.data = {}        
        self.plotting = Plotting(self)
        self.reporting = Reporting(self)
        self.log = Log(self)
        self.io = InputOutputHandling(self)
        self.solution_properties = SolutionProperties(self)
        
        # Several parts of the code wants to know these things,
        # so we keep them in a central place
        self.ndim = 0
        self.timestep = 0
        self.time = 0.0
        self.dt = 0.0
        self.restarted = False
        
        # These will be filled out when .setup() is configuring the Navier-Stokes
        # solver. Included here for documentation purposes only
        self.solver = None
        self.multi_phase_model = None
        self.mesh_morpher = None
        self.t_start = None
        self.probes = None
        
        # For timing the analysis
        self.prevtime = self.starttime = time.time()
        
    def setup(self):
        """
        Setup the simulation. This creates the .solver object as well as the mesh,
        boundary conditions, initial condition, function spaces, runtime
        post-processing probes, program and user defined hooks ... 
        """
        # The implementation is rather long, so it is in a separate file
        setup_simulation(self)
    
    def set_mesh(self, mesh, mesh_facet_regions=None):
        """
        Set the computational domain
        """
        self.data['mesh'] = mesh
        self.data['mesh_facet_regions'] = mesh_facet_regions
        self.ndim = mesh.topology().dim()
        self.update_mesh_data()
        self.log.info('Loaded mesh with %d cells' % mesh.num_cells())
    
    def update_mesh_data(self, connectivity_changed=True):
        """
        Some precomputed values must be calculated before the timestepping
        and updated every time the mesh changes
        """
        if connectivity_changed:
            init_connectivity(self)
        precompute_cell_data(self)
        precompute_facet_data(self)
        
        # Work around uflacs missing CellSize and CellVolume etc for isoparametric elements
        mesh = self.data['mesh']
        if mesh.ufl_coordinate_element().degree() > 1:
            # FIXME: this is only valid for uniform meshes!
            area = dolfin.assemble(1.0*dolfin.dx(mesh, degree=2))
            h = dolfin.Constant((area / mesh.num_cells())**(1.0 / mesh.topology().dim()))
        else:
            h = dolfin.CellSize(mesh)
        self.data['h'] = h
    
    def _at_start_of_timestep(self, timestep_number, t, dt):
        self.timestep = timestep_number
        self.time = t
        self.dt = dt
    
    def _at_end_of_timestep(self):
        # Report the time spent in this time step
        newtime = time.time()
        self.reporting.report_timestep_value('tstime', newtime-self.prevtime)
        self.reporting.report_timestep_value('tottime', newtime-self.starttime)
        self.prevtime = newtime
        
        # Report the solution properties
        if self.solution_properties.active:
            Co_max = self.solution_properties.courant_number().vector().max()
            Pe_max = self.solution_properties.peclet_number().vector().max()
            div_dS_f, div_dx_f = self.solution_properties.divergences()
            div_dS = div_dS_f.vector().max()
            div_dx = div_dx_f.vector().max()
            mass = self.solution_properties.total_mass()
            Ek, Ep = self.solution_properties.total_energy()
            self.reporting.report_timestep_value('Co', Co_max)
            self.reporting.report_timestep_value('Pe', Pe_max)
            self.reporting.report_timestep_value('div', div_dx + div_dS)
            self.reporting.report_timestep_value('mass', mass)
            self.reporting.report_timestep_value('Ek', Ek)
            self.reporting.report_timestep_value('Ep', Ep)
            
            if self.solution_properties.has_div_conv:
                # Convecting and convected velocities are separate
                div_conv_dS_f, div_conv_dx_f = self.solution_properties.divergences('u_conv')
                div_conv_dS = div_conv_dS_f.vector().max()
                div_conv_dx = div_conv_dx_f.vector().max()
                self.reporting.report_timestep_value('div_conv', div_conv_dx + div_conv_dS)
                
                # Difference between the convective and the convected velocity
                ucdiff = velocity_change(u1=self.data['up'],
                                         u2=self.data['up_conv'],
                                         ui_tmp=self.data['ui_tmp'])
                self.reporting.report_timestep_value('uconv_diff', ucdiff)
            
            if not numpy.isfinite(Co_max):
                ocellaris_error('Non finite Courant number',
                                'Found Co = %g' % Co_max)
        
        # Write fields to output file
        self.io.write_fields()
        
        # Write timestep report
        self.reporting.log_timestep_reports()
