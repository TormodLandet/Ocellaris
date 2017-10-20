"""
A simple elliptic solver for use in unit testing
"""
import dolfin
from dolfin import dot, grad, dx, dS, jump, avg
from ocellaris.solvers import Solver, register_solver
from ocellaris.solver_parts import define_penalty
from ocellaris.utils import linear_solver_from_input


BASE_INPUT = """
ocellaris:
    type: input
    version: 1.0

mesh:
    type: Rectangle
    Nx: 10
    Ny: 10

solver:
    type: PoissonDG
    phi:
        solver: cg
        preconditioner: hypre_amg
        parameters:
            relative_tolerance: 1.0e-10
            absolute_tolerance: 1.0e-15

output:
    solution_properties: off
    xdmf_write_interval: 0
    save_restart_file_at_end: off

# Dummy values below here are just to make Ocellaris run. The
# solution framework has some assumptions that all solvers
# are (quasi) time stepping flow solvers
time:
    dt: 1
physical_properties:
    nu0: 1
"""
 

@register_solver('PoissonDG')
class PoissonDGSolver(Solver):
    description = "Poisson equation solver using DG elements"
    
    @classmethod
    def create_function_spaces(cls, simulation):
        family = simulation.input.get_value('solver/function_space', 'DG', 'string')
        degree = simulation.input.get_value('solver/polynomial_degree', 1, 'int')
        mesh = simulation.data['mesh']
        simulation.data['Vphi'] = dolfin.FunctionSpace(mesh, family, degree)
    
    def __init__(self, simulation):
        """
        A discontinuous Galerkin Poisson solver for use in 
        the Ocellaris solution framework. Solves -∇⋅∇φ = f
        by use of the Symmetric Interior Penalty method
        """
        self.simulation = simulation
        self.setup_scalar_equation()
    
    def setup_scalar_equation(self):
        sim = self.simulation
        V = sim.data['Vphi']
        mesh = V.mesh()
        P = V.ufl_element().degree()
        
        # Source term
        source_cpp = sim.input.get_value('solver/source', '0', 'string')
        f = dolfin.Expression(source_cpp, degree=P)
        
        # Create the solution function 
        sim.data['phi'] = dolfin.Function(V)
        
        # DG elliptic penalty
        penalty = define_penalty(mesh, P, k_min=1.0, k_max=1.0)
        penalty_dS = dolfin.Constant(penalty)
        penalty_ds = dolfin.Constant(penalty*2)
        
        # Define weak form
        u, v = dolfin.TrialFunction(V), dolfin.TestFunction(V) 
        a = dot(grad(u), grad(v))*dx
        L = f*v*dx
        
        # Symmetric Interior Penalty method for -∇⋅∇φ
        n = dolfin.FacetNormal(mesh)
        a -= dot(n('+'), avg(grad(u)))*jump(v)*dS
        a -= dot(n('+'), avg(grad(v)))*jump(u)*dS
        
        # Symmetric Interior Penalty coercivity term
        a += penalty_dS*jump(u)*jump(v)*dS
        
        # Dirichlet boundary conditions
        dirichlet_bcs = sim.data['dirichlet_bcs'].get('phi', [])
        for dbc in dirichlet_bcs:
            bcval, dds = dbc.func(), dbc.ds()
            
            # SIPG for -∇⋅∇φ
            a -= dot(n, grad(u))*v*dds
            a -= dot(n, grad(v))*u*dds
            L -= dot(n, grad(v))*bcval*dds
            
            # Weak Dirichlet
            a += penalty_ds*u*v*dds
            L += penalty_ds*bcval*v*dds
        
        # Neumann boundary conditions
        neumann_bcs = sim.data['neumann_bcs'].get('phi', [])
        for nbc in neumann_bcs:
            L += nbc.func()*v*nbc.ds()
        
        # Does the system have a null-space?
        self.has_null_space = len(dirichlet_bcs) == 0
        
        self.form_lhs = a
        self.form_rhs = L
    
    def run(self):
        sim = self.simulation
        sim.hooks.simulation_started()
        sim.hooks.new_timestep(timestep_number=1, t=1.0, dt=1.0)
        
        # Assemble system
        A = dolfin.assemble(self.form_lhs)
        b = dolfin.assemble(self.form_rhs)
        
        # Create Krylov solver
        solver = linear_solver_from_input(sim, 'solver/phi', 'cg')
        solver.set_operator(A)
        
        # The function where the solution is stored
        phi = self.simulation.data['phi']
        
        # Remove null space if present (pure Neumann BCs)
        if self.has_null_space:
            null_vec = dolfin.Vector(phi.vector())
            phi.function_space().dofmap().set(null_vec, 1.0)
            null_vec *= 1.0/null_vec.norm("l2")
            
            null_space = dolfin.VectorSpaceBasis([null_vec])
            dolfin.as_backend_type(A).set_nullspace(null_space)
            null_space.orthogonalize(b)
        
        # Solve the linear system using the default solver (direct LU solver)
        solver.solve(phi.vector(), b)
        
        sim.hooks.end_timestep()
        sim.hooks.simulation_ended(success=True)
