import dolfin
from ocellaris.utils import timeit
from .io_impl import RestartFileIO, XDMFFileIO, DebugIO


# Default values, can be changed in the input file
XDMF_WRITE_INTERVAL = 0
HDF5_WRITE_INTERVAL = 0
SAVE_RESTART_AT_END = True


class InputOutputHandling():
    def __init__(self, simulation):
        """
        This class handles reading and writing the simulation state such as
        velocity and presure fields. Files for postprocessing (xdmf) are also
        handled here
        """
        self.simulation = sim = simulation
        self.ready = False
        self.persisted_python_data = {}
        
        self.restart = RestartFileIO(simulation, self.persisted_python_data)
        self.xdmf = XDMFFileIO(simulation)
        self.debug = DebugIO(simulation)
        
        sim.hooks.add_pre_simulation_hook(self._setup_io, 'Setup simulation IO')
        close = lambda _success: self._close_files()
        sim.hooks.add_post_simulation_hook(close, 'Save restart file and close files')
        
        # When exiting due to a signal kill/shutdown a savepoint file will be
        # written instead of an endpoint file, which is the default. This helps
        # with automatic restarts from checkpointing jobs. The only difference
        # is the name of the file
        self.last_savepoint_is_checkpoint = False
    
    def add_extra_output_function(self, function):
        """
        The output files (XDMF) normally only contain u, p and potentially rho or c. Other
        custom fields can be added  
        """
        self.simulation.log.info('    Adding extra output function %s' % function.name())
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
    
    def _setup_io(self):
        sim = self.simulation
        sim.log.info('Setting up simulation IO')
        
        # Make sure functions have nice names for output
        for name, description in (('p', 'Pressure'),
                                  ('p_hydrostatic', 'Hydrostatic pressure'),
                                  ('c', 'Colour function'),
                                  ('rho', 'Density'),
                                  ('u0', 'X-component of velocity'),
                                  ('u1', 'Y-component of velocity'),
                                  ('u2', 'Z-component of velocity'),
                                  ('boundary_marker', 'Domain boundary regions')):
            if not name in sim.data:
                continue
            func = sim.data[name]
            if hasattr(func, 'rename'):
                func.rename(name, description)
        
        # Dump initial state
        self.ready = True
        self.write_fields()
    
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
        elif sim.input.get_value('output/save_restart_file_at_end', SAVE_RESTART_AT_END, 'bool'):
            # Shutting down for good
            h5_file_name = sim.input.get_output_file_path('output/hdf5_file_name', '_endpoint_%08d.h5')
            h5_file_name = h5_file_name % sim.timestep
            self.write_restart_file(h5_file_name)
        
        self.xdmf.close()
    
    @timeit.named('IO write_fields')
    def write_fields(self):
        """
        Write fields to file after end of time step
        """
        sim = self.simulation
        
        # Allow these to change over time
        hdf5_write_interval = sim.input.get_value('output/hdf5_write_interval', HDF5_WRITE_INTERVAL, 'int')
        xdmf_write_interval = sim.input.get_value('output/xdmf_write_interval', XDMF_WRITE_INTERVAL, 'int')
        
        # No need to output just after restarting, this will overwrite the output from the previous simulation
        just_restarted = sim.restarted and sim.timestep_restart == 0
        
        if xdmf_write_interval > 0 and sim.timestep % xdmf_write_interval == 0 and not just_restarted:
            self.xdmf.write()
        
        if hdf5_write_interval > 0 and sim.timestep % hdf5_write_interval == 0 and not just_restarted:
            self.write_restart_file()
    
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
        except:
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
