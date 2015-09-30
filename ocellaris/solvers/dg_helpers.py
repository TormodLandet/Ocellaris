# encoding: utf8
import dolfin
from dolfin import MixedFunctionSpace, VectorFunctionSpace, FunctionSpace
from dolfin import FacetNormal, TrialFunction, TestFunctions, Function 
from dolfin import cells, dot, as_vector, dx, ds, dS, LocalSolver


def define_penalty(mesh, P, k_min, k_max, boost_factor=3, exponent=1):
    """
    Define the penalty parameter used in the Poisson equations
    
    Arguments:
        mesh: the mesh used in the simulation
        P: the polynomial degree of the unknown
        k_min: the minimum diffusion coefficient
        k_max: the maximum diffusion coefficient
        boost_factor: the penalty is multiplied by this factor
        exponent: set this to greater than 1 for superpenalisation
    """
    assert k_max >= k_min
    ndim = mesh.geometry().dim()
    
    # Calculate geometrical factor used in the penalty
    geom_fac = 0
    for cell in cells(mesh):
        vol = cell.volume()
        area = sum(cell.facet_area(i) for i in range(ndim + 1))
        gf = area/vol
        geom_fac = max(geom_fac, gf)
    
    penalty = boost_factor * k_max**2/k_min * (P + 1)*(P + ndim)/ndim * geom_fac**exponent
    return penalty


def bdm_projection(w, mesh):
    """
    Implement equation 4a and 4b in "Two new techniques for generating exactly
    incompressible approximate velocities" by Bernardo Cocburn (2009)
    
    For each element K in the mesh:
    
        <u⋅n, φ> = <û⋅n, φ>  ∀ ϕ ∈ P_{k}(F) for any face F ∈ ∂K
        (u, ϕ) = (w, ϕ)      ∀ φ ∈ P_{k-2}(K)^2
        (u, ϕ) = (w, ϕ)      ∀ φ ∈ {ϕ ∈ P_{k}(K)^2 : ∇⋅ϕ = 0 in K, ϕ⋅n = 0 on ∂K}  
        
    Here w is the input velocity function in DG2 space and û is the flux at
    each face (upwind value of the input velocity w). P_{x} is the space of
    polynomials of order k
    """
    ue = w[0].function_space().ufl_element()
    k = 2
    assert ue.family() == 'Discontinuous Lagrange'
    assert ue.degree() == k 
    assert w.ufl_shape == (2,)
    
    V = VectorFunctionSpace(mesh, 'DG', k)
    n = FacetNormal(mesh)
    
    # The mixed function space of the projection test functions
    V1 = FunctionSpace(mesh, 'DGT', k)
    V2 = VectorFunctionSpace(mesh, 'DG', k-2)
    V3 = FunctionSpace(mesh, 'Bubble', 3)
    W = MixedFunctionSpace([V1, V2, V3])
    v1, v2, v3b = TestFunctions(W)
    u = TrialFunction(V)
    
    # Calculate the upwind normal velocity, u_hat = w⋅n⁺, on all facets
    wn = dot(w, n)
    u_hat_ds = wn*n
    upwind = (wn + abs(wn))/2*n
    u_hat_dS = upwind('+') + upwind('-') 
    
    # Equation 1 - flux through the sides
    a = dot(u, n)*v1*ds
    L = dot(u_hat_ds, n)*v1*ds
    for R in '+-':
        a += dot(u(R), n(R))*v1(R)*dS
        L += dot(u_hat_dS, n(R))*v1(R)*dS 
    
    # Equation 2 - internal shape
    a += dot(u, v2)*dx
    L += dot(w, v2)*dx
    
    # Equation 3 - BDM Phi
    v3 = as_vector([v3b.dx(1), -v3b.dx(0)]) # Curl of [0, 0, v3b]
    a += dot(u, v3)*dx
    L += dot(w, v3)*dx
    
    # Find the projected velocity
    ls = LocalSolver(a, L)
    U = Function(V)
    ls.solve_local_rhs(U)
    
    a0 = dolfin.FunctionAssigner(w[0].function_space(), V.sub(0))
    a1 = dolfin.FunctionAssigner(w[1].function_space(), V.sub(1))
    
    U0, U1 = U.split()
    a0.assign(w[0], U0)
    a1.assign(w[1], U1)
