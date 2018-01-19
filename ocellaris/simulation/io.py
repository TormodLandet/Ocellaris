import os, re
import yaml
import numpy
import h5py
import dolfin
from ocellaris.utils import ocellaris_error, timeit


# Default values, can be changed in the input file
XDMF_WRITE_INTERVAL = 0
XDMF_FLUSH = True
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
        self.persisted_python_data = {}
        
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
        self.extra_xdmf_functions.append(function)
    
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
        self.xdmf_file = None
        
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
    
    def setup_xdmf(self):
        """
        Create XDMF file object
        """
        sim = self.simulation
        xdmf_flush =  sim.input.get_value('output/xdmf_flush', XDMF_FLUSH, 'bool')
        
        file_name = sim.input.get_output_file_path('output/xdmf_file_name', '.xdmf')
        file_name2 = os.path.splitext(file_name)[0] + '.h5'
        
        # Remove previous files
        if os.path.isfile(file_name) and sim.rank == 0:
            sim.log.info('    Removing existing XDMF file %s' % file_name)
            os.remove(file_name)
        if os.path.isfile(file_name2) and sim.rank == 0:
            sim.log.info('    Removing existing XDMF file %s' % file_name2)
            os.remove(file_name2)
        
        sim.log.info('    Creating XDMF file %s' % file_name)
        comm = sim.data['mesh'].mpi_comm()
        self.xdmf_file = dolfin.XDMFFile(comm, file_name)
        self.xdmf_file.parameters['flush_output'] = xdmf_flush
        self.xdmf_file.parameters['rewrite_function_mesh'] = False
        self.xdmf_file.parameters['functions_share_mesh'] = True
        self.xdmf_first_output = True
        
        def create_vec_func(V):
            "Create a vector function from the components"
            family = V.ufl_element().family()
            degree = V.ufl_element().degree()
            cd = sim.data['constrained_domain']
            V_vec = dolfin.VectorFunctionSpace(sim.data['mesh'], family, degree,
                                               constrained_domain=cd)
            vec_func = dolfin.Function(V_vec)
            assigner = dolfin.FunctionAssigner(V_vec, [V] * sim.ndim)
            return vec_func, assigner
        
        # XDMF cannot save functions given as "as_vector(list)" 
        self._vel_func, self._vel_func_assigner = create_vec_func(sim.data['Vu'])
        self._vel_func.rename('u', 'Velocity')
        if sim.mesh_morpher.active:
            self._mesh_vel_func, self._mesh_vel_func_assigner = create_vec_func(sim.data['Vmesh'])
            self._mesh_vel_func.rename('u_mesh', 'Velocity of the mesh')
    
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
        
        self.xdmf_file = None
    
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
            self.write_plot_file()
        
        if hdf5_write_interval > 0 and sim.timestep % hdf5_write_interval == 0 and not just_restarted:
            self.write_restart_file()
    
    def write_plot_file(self):
        """
        Write a file that can be used for visualization. The fluid fields will be automatically
        downgraded (interpolated) into something VTK can accept, typically linear CG elements.
        """
        with dolfin.Timer('Ocellaris save xdmf'):
            if self.xdmf_file is None:
                self.setup_xdmf()
            self._write_xdmf()
    
    def write_restart_file(self, h5_file_name=None):
        """
        Write a file that can be used to restart the simulation
        """
        with dolfin.Timer('Ocellaris save hdf5'):
            return self._write_hdf5(h5_file_name)
    
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
            self._read_hdf5(h5_file_name, read_input=True, read_results=False)
    
    def load_restart_file_results(self, h5_file_name):
        """
        Load the results stored on the given restart file
        """
        with dolfin.Timer('Ocellaris load hdf5'):
            self._read_hdf5(h5_file_name, read_input=False, read_results=True)
    
    def load_restart_file_functions(self, h5_file_name):
        """
        Load only the Functions stored on the given restart file
        Returns a dictionary of functions, does not affect the
        Simulation object itself (for switching meshes etc.)
        """
        with dolfin.Timer('Ocellaris load hdf5'):
            return self._read_hdf5_functions(h5_file_name)
    
    def _write_xdmf(self):
        """
        Write plot files for Paraview and similar applications
        """
        t = self.simulation.time
        
        if self.xdmf_first_output:
            bm = self.simulation.data['boundary_marker']
            self.xdmf_file.write(bm)
        
        # Write the fluid velocities
        self._vel_func_assigner.assign(self._vel_func, list(self.simulation.data['up']))
        self.xdmf_file.write(self._vel_func, t)
        
        # Write the mesh velocities (used in ALE calculations)
        if self.simulation.mesh_morpher.active:
            self._mesh_vel_func_assigner.assign(self._mesh_vel_func, list(self.simulation.data['u_mesh']))
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
        
        self.xdmf_first_output = False
    
    def _write_hdf5(self, h5_file_name=None):
        """
        Write fields to HDF5 file to support restarting the solver 
        """
        sim = self.simulation
        
        if h5_file_name is None:
            h5_file_name = sim.input.get_output_file_path('output/hdf5_file_name', '_savepoint_%08d.h5') 
            h5_file_name = h5_file_name % sim.timestep
        
        # Write dolfin objects using dolfin.HDF5File
        sim.log.info('Creating HDF5 restart file %s' % h5_file_name)
        comm = sim.data['mesh'].mpi_comm()
        with dolfin.HDF5File(comm, h5_file_name, 'w') as h5:
            # Write mesh
            h5.write(sim.data['mesh'], '/mesh')
            
            # Write mesh facet regions (from mesh generator)
            mfr = sim.data['mesh_facet_regions']
            if mfr is not None:
                h5.write(mfr, '/mesh_facet_regions')
            
            # Write functions, sorted to ensure deterministic output order
            funcnames = []
            skip = {'coupled', } # Skip these functions
            funcs = sorted(sim.data.items())
            for name, value in funcs:
                if isinstance(value, dolfin.Function) and name not in skip:
                    assert name not in funcnames
                    h5.write(value, '/%s' % name)
                    funcnames.append(name)
        
        # Only write metadata on root process
        # Important: no collective operations below this point!
        comm.barrier()
        if self.simulation.rank != 0:
            return h5_file_name
        
        # Write numpy objects and metadata using h5py.File
        with h5py.File(h5_file_name, 'r+') as hdf:
            # Metadata
            meta = hdf.create_group('ocellaris')
            meta.attrs['time'] = sim.time
            meta.attrs['iteration'] = sim.timestep
            meta.attrs['dt'] = sim.dt
            meta.attrs['restart_file_format'] = 2
            
            # Functions to save strings
            string_dt = h5py.special_dtype(vlen=str)
            def np_string(root, name, strdata):
                np_data = numpy.array(str(strdata).encode('utf8'), dtype=object)
                root.create_dataset(name, data=np_data, dtype=string_dt)
            def np_stringlist(root, name, strlist):
                np_list = numpy.array([str(s).encode('utf8') for s in strlist], dtype=object)
                root.create_dataset(name, data=np_list, dtype=string_dt)
            
            # List of names
            repnames = list(sim.reporting.timestep_xy_reports.keys())
            np_stringlist(meta, 'function_names', funcnames)
            np_stringlist(meta, 'report_names', repnames)
            
            # Save the current input and the full log file
            np_string(meta, 'input_file', sim.input)
            np_string(meta, 'full_log', sim.log.get_full_log())
            
            # Save reports
            reps = hdf.create_group('reports')
            reps['timesteps'] = numpy.array(sim.reporting.timesteps, dtype=float)
            for rep_name, values in sim.reporting.timestep_xy_reports.items():
                reps[rep_name] = numpy.array(values, dtype=float)
            
            # Save persistent data dictionaries
            pdd = hdf.create_group('ocellaris_data')
            i = 0
            for name, data in self.persisted_python_data.items():
                # Get stripped down data with only basic data types
                data2 = {}
                for k, v in data.items():
                    if is_basic_datatype(k) and is_basic_datatype(v):
                        data2[k] = v
                
                if not data2:
                    # no basic data to store, skip this dict
                    continue
                
                # Convert basic data type data to YAML format
                data3 = yaml.dump(data2)
                
                pdi = pdd.create_group('data_%02d' % i)
                np_string(pdi, 'name', name)
                np_string(pdi, 'data', data3)
                i += 1
        
        return h5_file_name
    
    def _read_hdf5_metadata(self, h5_file_name, function_details=False):
        """
        Read HDF5 restart file metadata
        """
        # Check file format and read metadata
        with h5py.File(h5_file_name, 'r') as hdf:
            if not 'ocellaris' in hdf:
                ocellaris_error('Error reading restart file',
                                'Restart file %r does not contain Ocellaris meta data'
                                % h5_file_name)
            
            meta = hdf['ocellaris']
            restart_file_version = meta.attrs['restart_file_format']
            if restart_file_version != 2:
                ocellaris_error('Error reading restart file',
                                'Restart file version is %d, this version of Ocellaris only ' %
                                restart_file_version + 'supports version 2')
            
            t = float(meta.attrs['time'])
            it = int(meta.attrs['iteration'])
            dt = float(meta.attrs['dt'])
            inpdata = meta['input_file'].value
            funcnames = list(meta['function_names'])
            
            # Read function signatures
            if function_details:
                signatures = {}
                for fname in funcnames:
                    signatures[fname] = hdf[fname].attrs['signature']
        
        if function_details:
            return funcnames, signatures
        else:
            return t, it, dt, inpdata, funcnames
    
    def _read_hdf5(self, h5_file_name, read_input=True, read_results=True):
        """
        Read an HDF5 restart file on the format written by _write_hdf5()
        """
        # Check file format and read metadata
        t, it, dt, inpdata, funcnames = self._read_hdf5_metadata(h5_file_name)
        
        sim = self.simulation
        h5 = dolfin.HDF5File(dolfin.MPI.comm_world, h5_file_name, 'r')
        
        # This flag is used in sim.setup() to to skip mesh creation
        # and may be used by user code etc
        sim.restarted = True
        
        if read_input:
            # Read the input file
            sim.input.read_yaml(yaml_string=inpdata)
            sim.input.file_name = h5_file_name
            sim.input.set_value('time/tstart', t)
            sim.input.set_value('time/dt', dt)
            
            # Read mesh data
            mesh = dolfin.Mesh()
            h5.read(mesh, '/mesh', False)
            if h5.has_dataset('/mesh_facet_regions'):
                mesh_facet_regions = dolfin.FacetFunction('size_t', mesh)
                h5.read(mesh_facet_regions, '/mesh_facet_regions')
            else:
                mesh_facet_regions = None
            sim.set_mesh(mesh, mesh_facet_regions)
        
        if read_results:
            sim.log.info('Reading fields from restart file %r' % h5_file_name)
            sim.timestep = it
            
            # Read result field functions
            for name in funcnames:
                sim.log.info('    Function %s' % name)
                h5.read(sim.data[name], '/%s' % name)
            
            h5.close() # Close dolfin.HDF5File
            
            # Read persisted data dictionaries with h5py
            with h5py.File(h5_file_name, 'r') as hdf:
                pdd = hdf.get('ocellaris_data', {})
                for pdi in pdd.values():
                    name = pdi['name'].value
                    data = pdi['data'].value
                    
                    # Parse YAML data and update the persistent data store
                    data2 = yaml.load(data)
                    data3 = self.get_persisted_dict(name)
                    data3.update(data2)
    
    def _read_hdf5_functions(self, h5_file_name):
        """
        Return a dictionary of Functions read from the HDF5 file
        Mixed or vector function spaces are not supported and
        None will be returned for these 
        """
        # Check file format and read metadata
        funcnames, signatures = self._read_hdf5_metadata(h5_file_name, function_details=True)
        
        # Read mesh data
        h5 = dolfin.HDF5File(dolfin.MPI.comm_world, h5_file_name, 'r')
        mesh = dolfin.Mesh()
        h5.read(mesh, '/mesh', False)
        
        def mk_func(name):
            signature = signatures[name].decode('utf8')
            
            # Parse strings like "FiniteElement('Discontinuous Lagrange', tetrahedron, 2)"
            pattern = r"FiniteElement\('(?P<family>[^']+)', \w+, (?P<degree>\d+)\)"
            m = re.match(pattern, signature)
            if not m:
                return None
            
            family = m.group('family')
            degree = int(m.group('degree'))
            
            self.simulation.log.info('        Found function %s of family %r '
                                     'and degree %d' % (name, family, degree))
            V = dolfin.FunctionSpace(mesh, family, degree)
            return dolfin.Function(V)
        
        # Read result field functions
        funcs = {}
        for name in funcnames:
            f = mk_func(name)
            if f is not None:
                h5.read(f, '/%s' % name)
            funcs[name] = f
        h5.close()
        return funcs
    
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
                for irow in range(value.size(0)):
                    indices, row_values = value.getrow(irow)
                    rows.append(len(indices) + rows[-1])
                    cols.extend(indices)
                    values.extend(row_values)
                data[key] = ('dolfin.PETScMatrix', (N, M), rows, cols, values)
            
            else:
                raise ValueError('Cannot save object of type %r' % type(value))
        
        import pickle
        with open(file_name, 'wb') as out:
            pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)
        self.simulation.log.info('Saved LA objects to %r (%r)' % (file_name, kwargs.keys()))
    
    
    def load_la_objects(self, file_name):
        """
        Load a dictionary of named PETScMatrix and PETScVector objects
        from a file
        """
        assert self.simulation.ncpu == 1, 'Not supported in parallel'
        
        import pickle
        with open(file_name, 'rb') as inp:
            data = pickle.load(inp)
        
        ret = {}
        for key, value in data.items():
            value_type = value[0]
            # Vectors
            if value_type == 'dolfin.PETScVector':
                arr = value[1]
                dolf_vec = dolfin.PETScVector(dolfin.MPI.comm_world, arr.size)
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


def is_basic_datatype(value):
    """
    We only save "basic" datatypes like ints, floats, strings and 
    lists, dictionaries or sets of these as pickles in restart files
    
    This function returns True if the datatype is such a basic data
    type  
    """
    if isinstance(value, (int, float, str, bytes, bool)):
        return True
    elif isinstance(value, (list, tuple, set)):
        return all(is_basic_datatype(v) for v in value)
    elif isinstance(value, dict):
        return all(is_basic_datatype(k) and is_basic_datatype(v) for k, v in value.items())
    return False
