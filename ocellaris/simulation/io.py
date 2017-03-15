import os
import numpy
import dolfin
from ocellaris.utils import ocellaris_error


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
        sim.hooks.add_pre_simulation_hook(self._setup_io, 'Setup simulation IO')
        close = lambda success: self._close_files()
        sim.hooks.add_post_simulation_hook(close, 'Save restart file and close files')
        self.extra_xdmf_functions = []
        
    def add_extra_output_function(self, function):
        """
        The output files (XDMF) normally only contain u, p and potentially rho or c. Other
        custom fields can be added  
        """
        self.simulation.log.info('    Adding extra output function %s' % function.name())
        self.extra_xdmf_functions.append(function)
    
    def _setup_io(self):
        sim = self.simulation
        sim.log.info('Setting up simulation IO')
        self.xdmf_write_interval = sim.input.get_value('output/xdmf_write_interval',
                                                       XDMF_WRITE_INTERVAL, 'int')
        self.hdf5_write_interval = sim.input.get_value('output/hdf5_write_interval',
                                                         HDF5_WRITE_INTERVAL, 'int')
        
        # Create XDMF file object
        create_vec_func = False
        if self.xdmf_write_interval > 0:
            create_vec_func = True
            file_name = sim.input.get_output_file_path('output/xdmf_file_name', '.xdmf')
            file_name2 = os.path.splitext(file_name)[0] + '.h5'
            
            # Remove previous files
            if os.path.isfile(file_name):
                sim.log.info('    Removing existing XDMF file %s' % file_name)
                os.remove(file_name)
            if os.path.isfile(file_name2):
                sim.log.info('    Removing existing XDMF file %s' % file_name2)
                os.remove(file_name2)
            
            sim.log.info('    Creating XDMF file %s' % file_name)
            self.xdmf_file = dolfin.XDMFFile(dolfin.mpi_comm_world(), file_name)
        
        def create_vec_func(V):
            "Create a vector function from the components"
            family = V.ufl_element().family()
            degree = V.ufl_element().degree()
            cd = sim.data['constrained_domain']
            V_vec = dolfin.VectorFunctionSpace(sim.data['mesh'], family, degree,
                                               constrained_domain=cd)
            vec_func = dolfin.Function(V_vec)
            
            # Create function assigners for the components
            assigners = [dolfin.FunctionAssigner(V_vec.sub(d), V) for d in range(sim.ndim)]
            
            return vec_func, assigners
        
        # Some output formats cannot save functions given as "as_vector(list)" 
        if create_vec_func:
            self._vel_func, self._vel_func_assigners = create_vec_func(sim.data['Vu'])
            self._vel_func.rename('u', 'Velocity')
        if sim.mesh_morpher.active and create_vec_func:
            self._mesh_vel_func, self._mesh_vel_func_assigners = create_vec_func(sim.data['Vmesh'])
            self._mesh_vel_func.rename('u_mesh', 'Velocity of the mesh')
        
        # Make sure functions have nice names for output
        for name, description in (('p', 'Pressure'),
                                  ('p_hydrostatic', 'Hydrostatic pressure'),
                                  ('c', 'Colour function'),
                                  ('rho', 'Density'),
                                  ('u0', 'X-component of velocity'),
                                  ('u1', 'Y-component of velocity'),
                                  ('u2', 'Z-component of velocity')):
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
        if sim.input.get_value('output/save_restart_file_at_end',
                               SAVE_RESTART_AT_END, 'bool'):
            h5_file_name = sim.input.get_output_file_path('output/hdf5_file_name', '_endpoint_%08d.h5') 
            h5_file_name = h5_file_name % sim.timestep
            self.write_restart_file(h5_file_name)
        
        if self.xdmf_write_interval > 0:
            del self.xdmf_file
    
    def write_fields(self):
        """
        Write fields to file after end of time step
        """
        sim = self.simulation
        
        if self.xdmf_write_interval > 0 and sim.timestep % self.xdmf_write_interval == 0:
            self.write_plot_file()
        
        if self.hdf5_write_interval > 0 and sim.timestep % self.hdf5_write_interval == 0:
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
        
        # Write the fluid velocities
        for d in range(self.simulation.ndim):
            ui = self.simulation.data['up%d' % d]
            self._vel_func_assigners[d].assign(self._vel_func.sub(d), ui)
        self.xdmf_file.write(self._vel_func, t)
        
        # Write the mesh velocities (used in ALE calculations)
        if self.simulation.mesh_morpher.active:
            for d in range(self.simulation.ndim):
                ui = self.simulation.data['u_mesh%d' % d]
                self._mesh_vel_func_assigners[d].assign(self._mesh_vel_func.sub(d), ui)
            self.xdmf_file.write(self._mesh_vel_func, t)
        
        # Write scalar functions
        for name in ('p', 'p_hydrostatic', 'c', 'rho'):
            if name in self.simulation.data:
                func = self.simulation.data[name]
                if isinstance(func, dolfin.Function): 
                    self.xdmf_file.write(func, t)
        
        # Write extra functions
        for func in self.extra_xdmf_functions:
            self.xdmf_file.write(func, t)
    
    def _write_hdf5(self, h5_file_name=None):
        """
        Write fields to HDF5 file to support restarting the solver 
        """
        sim = self.simulation
        
        if h5_file_name is None:
            h5_file_name = sim.input.get_output_file_path('output/hdf5_file_name', '_savepoint_%08d.h5') 
            h5_file_name = h5_file_name % sim.timestep
        
        # Create HDF5 file object
        sim.log.info('Creating HDF5 restart file %s' % h5_file_name)
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
        h5.attributes('/ocellaris')['functions'] = ','.join(funcnames)
        
        # Save the log taking into account that older HDF5 formats
        # have limits on attribute size 
        full_log = sim.log.get_full_log()
        N = len(full_log)
        M = 64*1000 # the HDF5 limit prior to 1.8.0
        i = 0
        while i*M < N:
            log_part = full_log[i*M:(i+1)*M]
            h5.attributes('/ocellaris')['full_log_%d' % i] = log_part
            i += 1
        
        # Save reports
        timesteps = numpy.array(sim.reporting.timesteps, dtype=float)
        h5.write(timesteps, '/reports/timesteps')
        for rep_name, values in sim.reporting.timestep_xy_reports.items():
            assert ',' not in rep_name
            values = numpy.array(values, dtype=float)
            h5.write(values, '/reports/%s' % rep_name)
        repnames = ','.join(sim.reporting.timestep_xy_reports)
        if repnames: # Size of a HDF5 string  value must be > 0
            h5.attributes('/reports')['report_names'] = repnames 
        
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
                mesh_facet_regions = dolfin.FacetFunction('size_t', mesh)
                h5.read(mesh_facet_regions, '/mesh_facet_regions')
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


    ###########################################################
    # Debugging routines
    # These routines are not used during normal runs
    
    def save_la_objects(self, file_name, **kwargs):
        """
        Save a dictionary of named PETScMatrix and PETScVector objects
        to a file
        """
        assert self.simulation.ncpu == 1, 'Not supported in parallel'
        
        data = {}
        for key, value in kwargs.items():
            # Vectors
            if isinstance(value, dolfin.PETScVector):
                data[key] = ('dolfin.PETScVector', value.get_local())
            # Matrices
            elif isinstance(value, dolfin.PETScMatrix):
                rows = [0]
                cols = []
                values = []
                N, M = value.size(0), value.size(1)
                for irow in xrange(value.size(0)):
                    indices, row_values = value.getrow(irow)
                    rows.append(len(indices) + rows[-1])
                    cols.extend(indices)
                    values.extend(row_values)
                data[key] = ('dolfin.PETScMatrix', (N, M), rows, cols, values)
            
            else:
                raise ValueError('Cannot save object of type %r' % type(value))
        
        import cPickle as pickle
        with open(file_name, 'wb') as out:
            pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)
        self.simulation.log.info('Saved LA objects to %r (%r)' % (file_name, kwargs.keys()))
    
    
    def load_la_objects(self, file_name):
        """
        Load a dictionary of named PETScMatrix and PETScVector objects
        from a file
        """
        assert self.simulation.ncpu == 1, 'Not supported in parallel'
        
        import cPickle as pickle
        with open(file_name, 'rb') as inp:
            data = pickle.load(inp)
        
        ret = {}
        for key, value in data.items():
            value_type = value[0]
            # Vectors
            if value_type == 'dolfin.PETScVector':
                arr = value[1]
                dolf_vec = dolfin.PETScVector(dolfin.mpi_comm_world(), arr.size)
                dolf_vec.set_local(arr)
                dolf_vec.apply('insert')
                ret[key] = dolf_vec
            # Matrices
            elif value_type == 'dolfin.PETScMatrix':
                shape, rows, cols, values = value[1:]                
                from petsc4py import PETSc
                mat = PETSc.Mat().createAIJ(size=shape, csr=(rows, cols, values))
                mat.assemble()
                dolf_mat = dolfin.PETScMatrix(mat)
                ret[key] = dolf_mat
            
            else:
                raise ValueError('Cannot save object of type %r' % value_type)

        self.simulation.log.info('Loaded LA objects from %r (%r)' % (file_name, data.keys()))
        return ret
