# encoding: utf8
import dolfin
from ocellaris.utils import timeit, ocellaris_interpolate


class MeshMorpher(object):
    def __init__(self, simulation):
        """
        Class to handle prescribed or runtime evaluated mesh morphing
        """
        self.simulation = simulation
        self.active = False
        
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
        assert self.active is False, 'Trying to setup mesh morphing twice in the same simulation'
        
        # Store previous cell volumes
        mesh = sim.data['mesh']
        Vcvol = dolfin.FunctionSpace(mesh, 'DG', 0) 
        sim.data['cvolp'] = dolfin.Function(Vcvol)
        
        # The function spaces for mesh velocities and displacements
        Vmesh = dolfin.FunctionSpace(mesh, 'CG', 1)
        Vmesh_vec = dolfin.VectorFunctionSpace(mesh, 'CG', 1)
        sim.data['Vmesh'] = Vmesh
        
        # Create mesh velocity functions
        u_mesh = []
        for d in range(sim.ndim):
            umi = dolfin.Function(Vmesh)
            sim.data['u_mesh%d' % d] = umi
            u_mesh.append(umi)
        u_mesh = dolfin.as_vector(u_mesh)
        sim.data['u_mesh'] = u_mesh
        
        # Divergence remover for mesh velocities
        self.divergence_remover = DivergenceRemover(u_mesh)
        
        # Create mesh displacement vector function
        self.displacement = dolfin.Function(Vmesh_vec)
        self.assigners = [dolfin.FunctionAssigner(Vmesh_vec.sub(d), Vmesh) for d in range(sim.ndim)]
        self.active = True
    
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
        Vmesh  = sim.data['Vmesh']
        
        for d, cpp_code in enumerate(self.cpp_codes):
            description = 'initial conditions for mesh vel %d' % d
            func = sim.data['u_mesh%d' % d]
            
            # Update the mesh velocity functions
            ocellaris_interpolate(sim, cpp_code, description, Vmesh, func)
        
        self.morph_mesh()
    
    def morph_mesh(self):
        """
        Move the mesh and update cached geometry information. It is assumed
        that the mesh velocities u_mesh0, u_mesh1 etc are already populated
        """
        sim = self.simulation
        
        # Remove any divergence from the velocity field
        self.divergence_remover.update_field()
        
        # Get the mesh displacement
        for d in range(sim.ndim):
            umi = sim.data['u_mesh%d' % d]
            self.assigners[d].assign(self.displacement.sub(d), umi)
        self.displacement.vector()[:] *= sim.dt
        
        # Save the cell volumes before morphing
        mesh = sim.data['mesh']
        cvolp = sim.data['cvolp']
        dofmap_cvol = cvolp.function_space().dofmap()
        for cell in dolfin.cells(mesh):
            dofs = dofmap_cvol.cell_dofs(cell.index())
            assert len(dofs) == 1
            cvolp.vector()[dofs[0]] = cell.volume()
        
        # Move the mesh according to the given displacements
        mesh.move(self.displacement)
        mesh.bounding_box_tree().build(mesh)
        sim.update_mesh_data(connectivity_changed=False)


class DivergenceRemover(object):
    def __init__(self, u_star):
        """
        Remove the divergence from a vector field by using Helmholtz decomposition.
        
            ∇⋅∇φ = ∇⋅u_star
            u_divfree = u_star - ∇φ
        
        The field is modified in place when .update_field() is called
        """
        self.u_star = u_star
        V = u_star[0].function_space()
        mesh = V.mesh()
        
        self.phi = dolfin.Function(V)
        self.uc = dolfin.Function(V)
        
        # Define weak form of Poisson problem
        from dolfin import grad, div, dot, dx, ds
        phi = dolfin.TrialFunction(V)
        q = dolfin.TestFunction(V)
        n = dolfin.FacetNormal(mesh)
        a = -dot(grad(phi), grad(q))*dx + dot(grad(phi), n)*q*ds
        L = div(u_star)*q*dx
        self.tensor_lhs = dolfin.assemble(a)
        self.form_rhs = L
        
        # Define velocity update equation
        self.tensor_lhs2 = dolfin.assemble(phi*q*dx)
        self.form_rhs2 = []
        for d in range(len(u_star)):
            self.form_rhs2.append(self.phi.dx(d)*q*dx)
    
    def update_field(self):
        """
        This function modifies u_star in place and removes the divergence
        """
        A = self.tensor_lhs
        b = dolfin.assemble(self.form_rhs)
        phi = self.phi
        u_star = self.u_star
        
        dolfin.solve(A, phi.vector(), b)
        
        uc = self.uc
        A2 = self.tensor_lhs2
        for d in range(len(u_star)):
            b2 = dolfin.assemble(self.form_rhs2[d])
            dolfin.solve(A2, uc.vector(), b2)
            u_star[d].vector().axpy(-1.0, uc.vector())
