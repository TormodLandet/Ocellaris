# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os
import dolfin
from ocellaris.utils import timeit
from .io_impl import RestartFileIO, XDMFFileIO, LegacyVTKIO, DebugIO


# Default values, can be changed in the input file
XDMF_WRITE_INTERVAL = 0
HDF5_WRITE_INTERVAL = 0
LVTK_WRITE_INTERVAL = 0
SAVE_RESTART_AT_END = True


class InputOutputHandling:
    def __init__(self, simulation):
        """
        This class handles reading and writing the simulation state such as
        velocity and presure fields. Files for postprocessing (xdmf) are also
        handled here
        """
        self.simulation = sim = simulation
        self.ready = False
        self.persisted_python_data = {}

        # Initialise the main plot file types
        self.restart = RestartFileIO(simulation, self.persisted_python_data)
        self.xdmf = XDMFFileIO(simulation)
        self.lvtk = LegacyVTKIO(simulation)
        self.debug = DebugIO(simulation)

        # Set up periodic output of plot files. Other parts of Ocellaris may
        # also add their plot writers to this list
        self._plotters = []
        self.add_plotter(
            self.xdmf.write, 'output/xdmf_write_interval', XDMF_WRITE_INTERVAL
        )
        self.add_plotter(
            self.lvtk.write, 'output/vtk_write_interval', LVTK_WRITE_INTERVAL
        )
        self.add_plotter(
            self._interval_write_restart,
            'output/hdf5_write_interval',
            HDF5_WRITE_INTERVAL,
        )

        def close(success):
            return self._close_files()  # @UnusedVariable - must be named success

        sim.hooks.add_post_simulation_hook(close, 'Save restart file and close files')

        # When exiting due to a signal kill/shutdown a savepoint file will be
        # written instead of an endpoint file, which is the default. This helps
        # with automatic restarts from checkpointing jobs. The only difference
        # is the name of the file
        self.last_savepoint_is_checkpoint = False

        # Last savepoint written. To be deleted after writing new save point
        # if only keeping the last save point file
        self.prev_savepoint_file_name = ''

    def setup(self):
        sim = self.simulation
        sim.log.info('Setting up simulation IO')

        # Make sure functions have nice names for output
        for name, description in (
            ('p', 'Pressure'),
            ('p_hydrostatic', 'Hydrostatic pressure'),
            ('c', 'Colour function'),
            ('rho', 'Density'),
            ('u0', 'X-component of velocity'),
            ('u1', 'Y-component of velocity'),
            ('u2', 'Z-component of velocity'),
            ('boundary_marker', 'Domain boundary regions'),
        ):
            if name not in sim.data:
                continue
            func = sim.data[name]
            if hasattr(func, 'rename'):
                func.rename(name, description)
        self.ready = True

    def add_plotter(self, func, interval_inp_key, default_interval):
        """
        Add a plotting function which produces IO output every timestep
        """
        self._plotters.append((func, interval_inp_key, default_interval))

    def add_extra_output_function(self, function):
        """
        The output files (XDMF) normally only contain u, p and potentially rho or c. Other
        custom fields can be added
        """
        self.simulation.log.info(
            '    Adding extra output function %s' % function.name()
        )
        self.xdmf.extra_functions.append(function)

    def get_persisted_dict(self, name):
        """
        Get dictionary that is persisted across program restarts by
        pickling the data when saving HDF5 restart files.

        Only basic data types in the dictionary are persisted. Such
        data types are ints, floats, strings, booleans and containers
        such as lists, dictionaries, tuples and sets of these basic
        data types. All other data can be stored in the returned
        dictionary, but will not be persisted
        """
        assert isinstance(name, str)
        return self.persisted_python_data.setdefault(name, {})

    def _close_files(self):
        """
        Save final restart file and close open files
        """
        if not self.ready:
            return

        sim = self.simulation
        if self.last_savepoint_is_checkpoint:
            # Shutting down, but ready to restart from checkpoint
            self.write_restart_file()
        elif sim.input.get_value(
            'output/save_restart_file_at_end', SAVE_RESTART_AT_END, 'bool'
        ):
            # Shutting down for good
            h5_file_name = sim.input.get_output_file_path(
                'output/hdf5_file_name', '_endpoint_%08d.h5'
            )
            h5_file_name = h5_file_name % sim.timestep
            self.write_restart_file(h5_file_name)

        self.xdmf.close()

    @timeit.named('IO write_fields')
    def write_fields(self):
        """
        Write fields to file after end of time step
        """
        sim = self.simulation

        # No need to output just after restarting, this will overwrite the
        # output from the previous simulation
        if sim.restarted and sim.timestep_restart == 0:
            return

        # Call the output functions at the right intervals
        for func, interval_inp_key, default_interval in self._plotters:
            # Check this every timestep, it might change
            write_interval = sim.input.get_value(
                interval_inp_key, default_interval, 'int'
            )

            # Write plot file this time step if it aligns with the interval
            if write_interval > 0 and sim.timestep % write_interval == 0:
                func()

    def _interval_write_restart(self):
        """
        Special treatment for restart files, to avoid filling up the disk we
        can delete the previous save point file after successfully writing
        a new save point file. This is only done for interval writes, not
        direct calls to write_restart_file ...
        """
        sim = self.simulation
        h5_file_name = self.write_restart_file()

        # Write was successfull (no exception) -> delete previous files
        if sim.input.get_value('output/hdf5_only_store_latest', False, 'bool'):
            if os.path.isfile(self.prev_savepoint_file_name) and sim.rank == 0:
                sim.log.info(
                    'Deleting previous save point file %r'
                    % self.prev_savepoint_file_name
                )
                os.unlink(self.prev_savepoint_file_name)
            self.prev_savepoint_file_name = h5_file_name

    def write_restart_file(self, h5_file_name=None):
        """
        Write a file that can be used to restart the simulation
        """
        with dolfin.Timer('Ocellaris save hdf5'):
            return self.restart.write(h5_file_name)

    def is_restart_file(self, file_name):
        """
        Is the given file an Ocellaris restart file
        """
        HDF5_SIGNATURE = b'\211HDF\r\n\032\n'
        try:
            # The HDF5 header is not guaranteed to be at offset 0, but for our
            # purposes this can be assumed as we do nothing special when writing
            # the HDF5 file (http://www.hdfgroup.org/HDF5/doc/H5.format.html).
            with open(file_name, 'rb') as inp:
                header = inp.read(8)
            return header == HDF5_SIGNATURE
        except Exception:
            return False

    def load_restart_file_input(self, h5_file_name):
        """
        Load the input used in the given restart file
        """
        with dolfin.Timer('Ocellaris load hdf5'):
            self.restart.read(h5_file_name, read_input=True, read_results=False)

    def load_restart_file_results(self, h5_file_name):
        """
        Load the results stored on the given restart file
        """
        with dolfin.Timer('Ocellaris load hdf5'):
            self.restart.read(h5_file_name, read_input=False, read_results=True)

    def load_restart_file_functions(self, h5_file_name):
        """
        Load only the Functions stored on the given restart file
        Returns a dictionary of functions, does not affect the
        Simulation object itself (for switching meshes etc.)
        """
        with dolfin.Timer('Ocellaris load hdf5'):
            return self.restart.read_functions(h5_file_name)
