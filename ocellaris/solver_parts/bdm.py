# encoding: utf8
import dolfin
from dolfin import FiniteElement, VectorElement, MixedElement, FunctionSpace, VectorFunctionSpace
from dolfin import FacetNormal, TrialFunction, TestFunctions, Function 
from dolfin import dot, as_vector, dx, ds, dS, LocalSolver


class VelocityBDMProjection():
    def __init__(self, simulation, w, incompressibility_flux_type='central', D12=None):
        """
        Implement equation 4a and 4b in "Two new techniques for generating exactly
        incompressible approximate velocities" by Bernardo Cocburn (2009)
        
        For each element K in the mesh:
        
            <u⋅n, φ> = <û⋅n, φ>  ∀ ϕ ∈ P_{k}(F) for any face F ∈ ∂K
            (u, ϕ) = (w, ϕ)      ∀ φ ∈ P_{k-2}(K)^2
            (u, ϕ) = (w, ϕ)      ∀ φ ∈ {ϕ ∈ P_{k}(K)^2 : ∇⋅ϕ = 0 in K, ϕ⋅n = 0 on ∂K}  
            
        Here w is the input velocity function in DG2 space and û is the flux at
        each face. P_{x} is the space of polynomials of order k

        The flux type can be 'central' or 'upwind'
        """
        ue = w[0].function_space().ufl_element()
        k = 2
        assert ue.family() == 'Discontinuous Lagrange'
        assert ue.degree() == k 
        assert w.ufl_shape == (2,)
        assert incompressibility_flux_type in ('central', 'upwind')
        
        dirichlet_bcs = simulation.data['dirichlet_bcs'].get('u', [])
        neumann_bcs = simulation.data['neumann_bcs'].get('u', [])
        
        mesh = w[0].function_space().mesh()
        V = VectorFunctionSpace(mesh, 'DG', k)
        n = FacetNormal(mesh)
        
        # The mixed function space of the projection test functions
        e1 = FiniteElement('DGT', mesh.ufl_cell(), k)
        e2 = VectorElement('DG', mesh.ufl_cell(), k-2)
        e3 = FiniteElement('Bubble', mesh.ufl_cell(), 3)
        em = MixedElement([e1, e2, e3])
        W = FunctionSpace(mesh, em)
        v1, v2, v3b = TestFunctions(W)
        u = TrialFunction(V)
        
        # The same fluxes that are used in the incompressibility equation
        if incompressibility_flux_type == 'central':    
            u_hat_dS = dolfin.avg(w)
        elif incompressibility_flux_type == 'upwind':
            w_nU = (dot(w, n) + abs(dot(w, n)))/2.0
            switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), 1.0, 0.0)
            u_hat_dS = switch*w('+') + (1 - switch)*w('-')
        
        if D12 is not None:
            u_hat_dS += dolfin.Constant([D12, D12])*dolfin.jump(w, n)
        
        # Equation 1 - flux through the sides
        a = dot(u, n)*v1*ds
        L = 0
        for R in '+-':
            a += dot(u(R), n(R))*v1(R)*dS
            L += dot(u_hat_dS, n(R))*v1(R)*dS
        for dbc in dirichlet_bcs:
            # FIXME: integrate only over inlet portions
            L += dot(dbc.func(), n)*v1*dbc.ds()
        for nbc in neumann_bcs:
            L += dot(w, n)*v1*nbc.ds()
         
        # Equation 2 - internal shape
        a += dot(u, v2)*dx
        L += dot(w, v2)*dx
        
        # Equation 3 - BDM Phi
        v3 = as_vector([v3b.dx(1), -v3b.dx(0)]) # Curl of [0, 0, v3b]
        a += dot(u, v3)*dx
        L += dot(w, v3)*dx
        
        # Pre-factorize matrices and store for usage in projection
        self.local_solver = LocalSolver(a, L)
        self.local_solver.factorize()
        self.temp_function = Function(V)
        self.w = w
        self.assigner0 = dolfin.FunctionAssigner(w[0].function_space(), V.sub(0))
        self.assigner1 = dolfin.FunctionAssigner(w[1].function_space(), V.sub(1))
    
    def run(self):
        """
        Perform the projection based on the current state of the Function w
        """
        # Find the projected velocity
        self.local_solver.solve_local_rhs(self.temp_function)
        
        # Assign to w
        U0, U1 = self.temp_function.split()
        self.assigner0.assign(self.w[0], U0)
        self.assigner1.assign(self.w[1], U1)
