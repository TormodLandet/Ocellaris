import numpy
import dolfin
from ocellaris.utils import ocellaris_error


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
            self.write_plot_file()
        
        if write_hdf5:
            self.write_restart_file()
            
    def write_plot_file(self):
        """
        Write a file that can be used for visualization. The fluid fields will be automatically
        downgraded (interpolated) into something VTK can accept, typically linear CG elements.
        """
        t = dolfin.Timer('Ocellaris save xdmf')
        self._write_xdmf()
        t.stop()
    
    def write_restart_file(self, h5_file_name=None):
        """
        Write a file that can be used to restart the simulation
        """
        t = dolfin.Timer('Ocellaris save hdf5')
        self._write_hdf5(h5_file_name)
        t.stop()
    
    def is_restart_file(self, file_name):
        """
        Is the given file an Ocellaris restart file
        """
        HDF5_SIGNATURE = '\211HDF\r\n\032\n'
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
        t = dolfin.Timer('Ocellaris load hdf5')
        self._read_hdf5(h5_file_name, read_input=True, read_results=False)
        t.stop()
        
    def load_restart_file_results(self, h5_file_name):
        """
        Load the results stored on the given restart file
        """
        t = dolfin.Timer('Ocellaris load hdf5')
        self._read_hdf5(h5_file_name, read_input=False, read_results=True)
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
    
    def _write_hdf5(self, h5_file_name=None):
        """
        Write fields to HDF5 file to support restarting the solver 
        """
        sim = self.simulation
        
        if h5_file_name is None:
            h5_file_name = sim.input.get_output_file_path('output/hdf5_file_name', '_savepoint_%r.h5') 
            h5_file_name = h5_file_name % sim.time
        
        # Create HDF5 file object
        h5 = dolfin.HDF5File(dolfin.mpi_comm_world(), h5_file_name, 'w')
        
        # Skip these functions
        skip = {'coupled', }
        
        # Write mesh
        h5.write(sim.data['mesh'], '/mesh')
        if sim.data['mesh_facet_regions'] is not None:
            h5.write(sim.data['mesh_facet_regions'], '/mesh_facet_regions')
            
        # Write functions
        funcnames = []
        for name, value in sim.data.items():
            if isinstance(value, dolfin.Function) and name not in skip:
                h5.write(value, '/%s' % name)
                
                # Save function names in a separate HDF attribute due to inability to 
                # list existing HDF groups when using the dolfin HDF5Function wrapper 
                assert ',' not in name
                funcnames.append(name) 
        
        # Metadata
        tinfo = numpy.array([sim.time, sim.timestep, sim.dt])
        h5.write(tinfo, '/ocellaris/time_info')
        h5.attributes('/ocellaris')['time'] = sim.time
        h5.attributes('/ocellaris')['iteration'] = sim.timestep
        h5.attributes('/ocellaris')['restart_file_format'] = 1
        h5.attributes('/ocellaris')['input_file'] = str(sim.input)
        h5.attributes('/ocellaris')['full_log'] = sim.log.get_full_log()
        h5.attributes('/ocellaris')['functions'] = ','.join(funcnames)
        
        h5.close()
        
    def _read_hdf5(self, h5_file_name, read_input=True, read_results=True):
        """
        Read an HDF5 restart file on the format written by _write_hdf5()
        """
        # Check for valid restart file
        h5 = dolfin.HDF5File(dolfin.mpi_comm_world(), h5_file_name, 'r')
        if not h5.has_dataset('ocellaris'):
            ocellaris_error('Error reading restart file',
                            'Restart file %r does not contain Ocellaris meta data'
                            % h5_file_name)
        restart_file_version = h5.attributes('/ocellaris')['restart_file_format']
        if restart_file_version != 1:
            ocellaris_error('Error reading restart file',
                            'Restart file version is %d, this version of Ocellaris only ' %
                            restart_file_version + 'supports version 1')
        
        # Read ocellaris metadata from h5 restart file
        t = h5.attributes('/ocellaris')['time']
        it = h5.attributes('/ocellaris')['iteration']
        inpdata = h5.attributes('/ocellaris')['input_file']
        funcnames = h5.attributes('/ocellaris')['functions'].split(',')
        
        sim = self.simulation
        
        if read_input:
            # Read the input file
            sim.input.read_yaml(yaml_string=inpdata)
            sim.input.set_value('time/tstart', t)
            
            # Read mesh data
            mesh = dolfin.Mesh()
            h5.read(mesh, '/mesh', False)
            if h5.has_dataset('/mesh_facet_regions'):
                mesh_facet_regions = dolfin.FacetFunction()
                h5.read(mesh_facet_regions, '/mesh_facet_regions', False)
            else:
                mesh_facet_regions = None
            sim.set_mesh(mesh, mesh_facet_regions)
            
            # This flag is used in sim.setup() to to skip mesh creation
            sim.restarted = True
        
        if read_results:
            sim.log.info('Reading fields from restart file %r' % h5_file_name)
            sim.timestep = it
            
            # Read result field functions
            for name in funcnames:
                sim.log.info('    Function %s' % name)
                h5.read(sim.data[name], '/%s' % name)
