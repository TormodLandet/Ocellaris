import dolfin

# Default values, can be changed in the input file
XDMF_WRITE_INTERVAL = 0
HDF5_WRITE_INTERVAL = 0

class InputOutputHandling():
    def __init__(self, simulation):
        """
        This class handles reading and writing the simulation state such as
        velocity and presure fields. Files for postprocessing (xdmf) are also
        handled here
        """
        self.simulation = sim = simulation
        sim.hooks.add_pre_simulation_hook(self._setup_io, 'Setup simulation IO')
    
    def _setup_io(self):
        sim = self.simulation
        self.xdmf_write_interval = sim.input.get_value('output/xdmf_write_interval',
                                                       XDMF_WRITE_INTERVAL, 'int')
        self.hdf5_write_interval = sim.input.get_value('output/hdf5_write_interval',
                                                         HDF5_WRITE_INTERVAL, 'int')
        
        # Create XDMF file object
        create_vec_func = False
        if self.xdmf_write_interval > 0:
            create_vec_func = True
            file_name = sim.input.get_output_file_path('output/xdmf_file_name', '.xdmf')
            self.xdmf_file = dolfin.XDMFFile(dolfin.mpi_comm_world(), file_name)
        
        # Create a vector function from the components
        if create_vec_func:
            Vu = sim.data['Vu']
            family = Vu.ufl_element().family()
            degree = Vu.ufl_element().degree()
            cd = sim.data['constrained_domain']
            Vu_vec = dolfin.VectorFunctionSpace(sim.data['mesh'], family, degree,
                                            constrained_domain=cd)
            self._vel_func = dolfin.Function(Vu_vec)
            self._vel_func.rename('u', 'Velocity')
            
            # Create function assigners for the components
            self._vel_func_assigners = [dolfin.FunctionAssigner(Vu_vec.sub(d), Vu) for d in range(sim.ndim)]
        
        # Make sure functions have nice names for output
        for name, description in (('p', 'Pressure'),
                                  ('p_hydrostatic', 'Hydrostatic pressure'),
                                  ('c', 'Colour function'),
                                  ('u0', 'X-component of velocity'),
                                  ('u1', 'Y-component of velocity'),
                                  ('u2', 'Z-component of velocity')):
            if name in sim.data:
                sim.data[name].rename(name, description)
    
    def write_fields(self):
        """
        Write fields to file after end of time step
        """
        sim = self.simulation
        
        # Should we write to file at all this timestep?
        write_xdmf = self.xdmf_write_interval > 0 and sim.timestep % self.xdmf_write_interval == 0
        write_hdf5 = self.hdf5_write_interval > 0 and sim.timestep % self.hdf5_write_interval == 0
        if not (write_xdmf or write_hdf5):
            return
        
        # Copy the velocity components into a vector function
        for d in range(sim.ndim):
            ui = sim.data['u%d' % d]
            self._vel_func_assigners[d].assign(self._vel_func.sub(d), ui)
        
        if write_xdmf:
            t = dolfin.Timer('Ocellaris save xdmf')
            self._write_xdmf()
            t.stop()
        
        if write_hdf5:
            self.write_restart_file()
    
    def write_restart_file(self):
        """
        Write a file that can be used to restart the simulation
        """
        t = dolfin.Timer('Ocellaris save hdf5')
        self._write_hdf5()
        t.stop()
    
    def _write_xdmf(self):
        """
        Write plot files for Paraview and similar applications
        """
        t = self.simulation.time
        self.xdmf_file << (self._vel_func, t)
        for name in ('p', 'p_hydrostatic', 'c'):
            if name in self.simulation.data:
                func = self.simulation.data[name] 
                self.xdmf_file << (func, t)
    
    def _write_hdf5(self):
        """
        Write fields to HDF5 file to support restarting the solver 
        """
        raise NotImplementedError('Untested code ...')
        
        sim = self.simulation
        
        # Create HDF5 file object
        h5_file_name = sim.input.get_output_file_path('output/hdf5_file_name', '_savepoint_%r.h5') 
        h5_file_name = h5_file_name % sim.time
        h5 = dolfin.HDF5File(dolfin.mpi_comm_world(), h5_file_name, 'w')
        
        h5.write(sim.data['mesh'], '/mesh')
        
        for name, value in sim.data.items():
            if isinstance(value, dolfin.Function):
                h5.write(value, '/%s' % name)
        
        # TODO: Save timestep, etc on the HDF5 file. Maybe also input 
        #       dictionary to restart directly from the h5 file?
        #       This would probably require h5py
        
        h5.close()
