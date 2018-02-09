import os
import dolfin


# Default values, can be changed in the input file
XDMF_FLUSH = True


class XDMFFileIO():
    def __init__(self, simulation):
        """
        Xdmf output using dolfin.XDMFFile
        """
        self.simulation = simulation
        self.extra_functions = []
        self.xdmf_file = None
    
    def close(self):
        """
        Close any open xdmf file
        """
        if self.xdmf_file is not None:
            self.xdmf_file.close()
        self.xdmf_file = None
    
    def write(self):
        """
        Write a file that can be used for visualization. The fluid fields will
        be automatically downgraded (interpolated) into something that
        dolfin.XDMFFile can write, typically linear CG elements.
        """
        with dolfin.Timer('Ocellaris save xdmf'):
            if self.xdmf_file is None:
                self._setup_xdmf()
            self._write_xdmf()
    
    def _setup_xdmf(self):
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
        for func in self.extra_functions:
            self.xdmf_file.write(func, t)
        
        self.xdmf_first_output = False
