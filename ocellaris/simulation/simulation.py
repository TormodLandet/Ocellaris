# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import time
import dolfin
from ocellaris.utils import timeit, shift_fields, verify_field_variable_definition
from ocellaris.utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data
from .hooks import Hooks
from .input import Input
from .reporting import Reporting
from .log import Log
from .io import InputOutputHandling
from .solution_properties import SolutionProperties
from .setup import setup_simulation


# Flush log and other output files at regular intervals, but not every
# timestep in case there are a lot of them per second
FLUSH_INTERVAL = 5  # seconds


class Simulation(object):
    def __init__(self):
        """
        Represents one Ocellaris simulation. The Simulation class
        connects the input file, geometry, mesh and more with the
        solver, results IO and reporting tools
        """
        # COMM_WORLD rank and size (may not match the mesh.mpi_comm())
        self.ncpu = dolfin.MPI.size(dolfin.MPI.comm_world)
        self.rank = dolfin.MPI.rank(dolfin.MPI.comm_world)

        self.data = {}  # Unknowns, coefficients and solutions
        self.fields = {}  # Known fields (incoming waves etc)

        self.hooks = Hooks(self)
        self.input = Input(self)
        self.reporting = Reporting(self)
        self.log = Log(self)
        self.io = InputOutputHandling(self)
        self.solution_properties = SolutionProperties(self)

        # Several parts of the code wants to know these things,
        # so we keep them in a central place
        self.ndim = 0
        self.timestep = 0  # Number of timesteps since beginning of sim
        self.timestep_restart = 0  # Number of timesteps since last restart
        self.time = 0.0
        self.dt = 0.0
        self.dt_prev = 0.0  # Some solvers need to know when the timestep changes
        self.restarted = False  # Starting from a restart file or from inp
        self.ndofs = 0

        # These will be filled out when .setup() is configuring the Navier-Stokes
        # solver. Included here for documentation purposes only
        self.solver = None
        self.multi_phase_model = None
        self.mesh_morpher = None
        self.t_start = None
        self.probes = None
        self.iso_surface_locators = {}

        # For timing the analysis and flushing the log at intervals
        self.prevtime = self.starttime = time.time()
        self.prevflush = 0

    @timeit.named('setup simulation')
    def setup(self):
        """
        Setup the simulation. This creates the .solver object as well as the mesh,
        boundary conditions, initial condition, function spaces, runtime
        post-processing probes, program and user defined hooks ...
        """
        self.flush_interval = self.input.get_value('output/flush_interval', FLUSH_INTERVAL, 'float')
        setup_simulation(self)

    def set_mesh(self, mesh, mesh_facet_regions=None):
        """
        Set the computational domain
        """
        self.data['mesh'] = mesh
        self.data['mesh_facet_regions'] = mesh_facet_regions
        self.ndim = mesh.topology().dim()
        assert self.ndim == mesh.geometry().dim()
        self.update_mesh_data()

        num_cells_local = mesh.topology().ghost_offset(self.ndim)
        num_cells_tot = dolfin.MPI.sum(mesh.mpi_comm(), float(num_cells_local))
        num_cells_min = dolfin.MPI.min(mesh.mpi_comm(), float(num_cells_local))
        num_cells_max = dolfin.MPI.max(mesh.mpi_comm(), float(num_cells_local))
        n_proc_mesh = dolfin.MPI.size(mesh.mpi_comm())
        self.log.info('Loaded mesh with %d cells' % num_cells_tot)
        self.log.info('    Distributed over %d MPI processes' % n_proc_mesh)
        self.log.info('    Least loaded process has %d cells' % num_cells_min)
        self.log.info('    Most loaded process has %d cells' % num_cells_max)

    def update_mesh_data(self, connectivity_changed=True):
        """
        Some precomputed values must be calculated before the timestepping
        and updated every time the mesh changes
        """
        if connectivity_changed:
            init_connectivity(self)
        precompute_cell_data(self)
        precompute_facet_data(self)

        # Work around missing consensus on what CellDiameter is for bendy cells
        mesh = self.data['mesh']
        if mesh.ufl_coordinate_element().degree() > 1:
            h = dolfin.CellVolume(mesh) ** (1 / mesh.topology().dim())
        else:
            h = dolfin.CellDiameter(mesh)
        self.data['h'] = h

    def get_data(self, name):
        """
        Return a solver variable if one exists with the given name (p, u0,
        mesh, ...) or a known field function if no solver variable with
        the given name exists.
        
        Known field functions must be specified (as always) with forward
        slash separated field name and function name, e.g., "waves/c"
        """
        if name in self.data:
            return self.data[name]
        else:
            return verify_field_variable_definition(self, name, 'Simulation.get_data()')

    def _at_start_of_simulation(self):
        """
        Runs after all the pre_simulation hooks, before staring the solver
        """
        # Give reasonable starting guesses for the solvers and something
        # sensible to use when computing the initial solution properties
        if 'up0' in self.data:
            shift_fields(self, ['up%d', 'u%d'])

        # Show solution properties without adding to the timestep reports
        # This should be before write_fields in case plot_divergences is on
        self.solution_properties.report(create_report=False)

        # Setup IO and dump initial fields
        self.io.setup()
        self.io.write_fields()

        self.flush(force=True)

    def _at_end_of_simulation(self, success):
        """
        Runs after the solver and after any post_simulation hooks
        """
        self.flush(force=True)

    def _at_start_of_timestep(self, timestep_number, t, dt):
        """
        Runs before any pre_timestep hooks, at the very beginning of the each
        time step in the solver
        """
        self.timestep = timestep_number
        self.timestep_restart += 1
        self.time = t
        self.dt_prev = self.dt
        self.dt = dt

        if 'dt' in self.data:
            self.data['dt'].assign(dt)

    @timeit.named('simulation at_end_of_timestep')
    def _at_end_of_timestep(self):
        """
        Runs after each time step in the solver, after all post_timestep hooks
        """
        # Report the time spent in this time step
        newtime = time.time()
        self.reporting.report_timestep_value('tstime', newtime - self.prevtime)
        self.reporting.report_timestep_value('tottime', newtime - self.starttime)
        self.prevtime = newtime

        # Write timestep report
        self.solution_properties.report()
        self.reporting.log_timestep_reports()

        # Write fields to output file
        self.io.write_fields()

        self.flush()

    def flush(self, force=False):
        """
        Flush output files if an appropriate amount of time has passed. This
        ensures that flush can be called after important output without slowing
        down the solver too much with disk IO in case of many calls to flush in
        quick succession
        """
        now = time.time()
        if force or now - self.prevflush > self.flush_interval:
            self.hooks.run_custom_hook('flush')
            self.prevflush = now
