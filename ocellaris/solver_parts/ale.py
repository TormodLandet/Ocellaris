import dolfin
from ocellaris.utils import timeit, ocellaris_project

class MeshMorpher(object):
    def __init__(self, simulation):
        """
        Class to handle prescribed or runtime evaluated mesh morphing
        """
        self.simulation = simulation
        
        # The user can give a mesh velocity function to simulate a piston or similar 
        prescribed_velocity_input = simulation.input.get_value('mesh/prescribed_velocity', None)
        if prescribed_velocity_input is not None:
            self.setup_prescribed_velocity(prescribed_velocity_input)
            
        # Fields that will be interpolated to the new location
        self.previous_fields = ['up0', 'up1', 'up2',
                                'upp0', 'upp1', 'upp2',
                                'uppp0', 'uppp1', 'uppp2']
    
    def setup(self):
        """
        Create mesh velocity and deformation functions
        """
        sim = self.simulation
        assert 'u_mesh' not in sim.data
        
        # The function spaces for mesh velocities and displacements
        mesh = sim.data['mesh']
        Vu = sim.data['Vu']
        Vmesh = dolfin.FunctionSpace(mesh, 'CG', 1)
        Vmesh_vec = dolfin.VectorFunctionSpace(mesh, 'CG', 1)
        
        # Create mesh velocity functions
        u_mesh = []
        for d in range(sim.ndim):
            umi = dolfin.Function(Vu)
            sim.data['u_mesh%d' % d] = umi
            u_mesh.append(umi)
        sim.data['u_mesh'] = dolfin.as_vector(u_mesh)
        
        # Create mesh displacement vector function
        self.displacement = dolfin.Function(Vmesh_vec)
        self.assigners = [dolfin.FunctionAssigner(Vmesh_vec.sub(d), Vmesh) for d in range(sim.ndim)]
        self.displacement_component = dolfin.Function(Vmesh)
    
    def setup_prescribed_velocity(self, input_dict):
        """
        The mesh nodes will be updated onece every time step 
        """
        # Read and verify input
        assert input_dict['type'] == 'CppCodedValue' 
        self.cpp_codes = input_dict['cpp_code']
        assert isinstance(self.cpp_codes, list) and len(self.cpp_codes) == self.simulation.ndim 
        
        # Setup mesh morphing every time step
        self.setup()
        self.simulation.hooks.add_pre_timestep_hook(lambda it, t, dt: self.update_prescribed_mesh_velocity(),
                                                    'MeshMorpher - update prescribed mesh velocity')
    
    @timeit
    def update_prescribed_mesh_velocity(self):
        """
        Move the mesh according to prescribed velocities
        """
        sim = self.simulation
        Vu  = sim.data['Vu']
        
        for d, cpp_code in enumerate(self.cpp_codes):
            description = 'initial conditions for mesh vel %d' % d
            func = sim.data['u_mesh%d' % d]
            
            # Update the mesh velocity functions
            ocellaris_project(sim, cpp_code, description, Vu, func)
        
        self.morph_mesh()
    
    def morph_mesh(self):
        """
        Move the mesh and update cached geometry information. It is assumed
        that the mesh velocities u_mesh0, u_mesh1 etc are already populated
        """
        sim = self.simulation
        dc = self.displacement_component
        Vmesh = dc.function_space()
        
        # Get the mesh displacement
        for d in range(sim.ndim):
            umi = sim.data['u_mesh%d' % d]
            dolfin.project(umi, V=Vmesh, function=self.displacement_component)
            self.displacement_component.vector()[:] *= sim.dt
            self.assigners[d].assign(self.displacement.sub(d), self.displacement_component)
            
        # Save the old mesh for interpolation
        mesh = sim.data['mesh']
        old_mesh = dolfin.Mesh(mesh)
        
        # Move the mesh according to the given displacements
        mesh.move(self.displacement)
        sim.update_mesh_data(connectivity_changed=False)
        
        # Interpolate previous fields to the new mesh location
        for name in self.previous_fields:
            if not name in sim.data:
                continue
            
            # The field has the new mesh, but the values have not been updated
            field_new = sim.data[name]
            V_new = field_new.function_space()
            
            # Make a field with the old mesh and use the same dof values
            family = V_new.ufl_element().family()
            degree = V_new.ufl_element().degree()
            V_old = dolfin.FunctionSpace(old_mesh, family, degree)
            field_old = dolfin.Function(V_old)
            field_old.assign(field_new)
            
            # Interpolate to the new locations
            interpolated = dolfin.interpolate(field_old, V_new)
            field_new.assign(interpolated)
