from dolfin import dot, grad, dx, assemble, Timer, TrialFunction, TestFunction
from dolfin import PETScKrylovSolver, sqrt, VectorSpaceBasis, as_backend_type, PETScOptions


class HydrostaticPressure(object):
    def __init__(self, rho, g, ph):
        """
        Calculate the hydrostatic pressure

        The gravity vector g *must* be parallel to one of the axes
        """
        self.active = True
        
        # Check if there is any gravity
        if all(gi == 0 for gi in g.values()):
            self.active = False
            return
        
        # Define the weak form
        Vp = ph.function_space()
        p = TrialFunction(Vp)
        q = TestFunction(Vp)
        a = dot(grad(p), grad(q))*dx
        L = rho*dot(g, grad(q))*dx
        
        self.func = ph
        self.tensor_lhs = assemble(a)
        self.form_rhs = L
        self.null_space = None
        self.solver = None
    
    def update(self):
        if not self.active:
            self.func.zero()
            return
        
        t = Timer('Ocellaris update hydrostatic pressure')
        
        A = self.tensor_lhs
        b = assemble(self.form_rhs)
        
        if self.solver is None:
            self.null_space, self.solver = setup_solver(A, b)
        
        self.null_space.orthogonalize(b)
        self.solver.solve(A, self.func.vector(), b)
        
        t.stop()


def setup_solver(A, b, tol=1e-15, verbose=False):
    """
    This is taken from the Dolfin singular Poisson demo
    with a few modifications to avoid stomping on other
    solvers KSPOptions (by setting a unique prefix)
    """
    PREFIX = 'ocellaris_hydrostatic_pressure_'
    solver = PETScKrylovSolver()
    solver.set_options_prefix(PREFIX)
    
    # Create vector that spans the null space, and normalize
    null_space_vector = b.copy()
    null_space_vector[:] = sqrt(1.0/null_space_vector.size())
    
    # Create null space basis object and attach to PETSc matrix
    null_space = VectorSpaceBasis([null_space_vector])
    as_backend_type(A).set_nullspace(null_space)
    
    # Set PETSc solve type (conjugate gradient) and preconditioner
    # (algebraic multigrid)
    PETScOptions.set(PREFIX + "ksp_type", "cg")
    PETScOptions.set(PREFIX + "pc_type", "gamg")
    
    # Since we have a singular problem, use SVD solver on the multigrid
    # 'coarse grid'
    PETScOptions.set(PREFIX + "mg_coarse_ksp_type", "preonly")
    PETScOptions.set(PREFIX + "mg_coarse_pc_type", "svd")
    
    # Set the solver tolerance and allow starting from previous solution
    PETScOptions.set(PREFIX + "ksp_rtol", tol)
    PETScOptions.set(PREFIX + "ksp_initial_guess_nonzero", True)
    
    if verbose:
        # Print PETSc solver configuration
        PETScOptions.set(PREFIX + "ksp_view")
        PETScOptions.set(PREFIX + "ksp_monitor")
    
    # Create Krylov solver and set operator
    solver.set_operator(A)
    
    # Set PETSc options on the solver
    solver.set_from_options()
    
    return null_space, solver
