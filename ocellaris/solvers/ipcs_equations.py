# encoding: utf8
import dolfin
from dolfin import dot, nabla_grad, avg, jump, dx, dS

def define_advection_problem(u, v, up, upp, u_conv1, f, n, beta, time_coeffs, dt, dirichlet_bcs):
    """
    Define the advection problem
    
     d/dt(u) + u_conv1 ⋅ grad(u) = f
     
    Returns the bilinear and linear forms
    """
    family = u.element().family()
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        c1, c2, c3 = time_coeffs 
        eq = (c1*u + c2*up + c3*upp)/dt*v*dx + dot(u_conv1, nabla_grad(u))*v*dx - f*v*dx
    
    elif family == 'Discontinuous Lagrange':
        # Upstream and downstream normal velocities
        flux_nU = u*(dot(u_conv1, n) + abs(dot(u_conv1, n)))/2
        flux_nD = u*(dot(u_conv1, n) - abs(dot(u_conv1, n)))/2
        
        # Define the blended flux
        # The blending factor beta is not DG, so beta('+') == beta('-')
        b = beta('+')
        flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
        
        # Equation to solve
        c1, c2, c3 = time_coeffs 
        eq = (c1*u + c2*up + c3*upp)/dt*v*dx \
             - u*dot(u_conv1, nabla_grad(v))*dx \
             + flux*jump(v)*dS  - f*v*dx
        
        # Enforce Dirichlet BCs weakly
        for dbc in dirichlet_bcs:
            eq += dot(u_conv1, n)*v*(u - dbc.func())*dbc.ds()
        
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)    
    return a, L

def define_poisson_problem(u, v, k, f, n, penalty, dirichlet_bcs, neumann_bcs):
    """
    Define the Poisson problem for u in f.space V
    
        - div(k*grad(u)) = f
    
    Note the minus in front of the first term!
    
    Arguments:
        u: a TrialFunction
        v: a TestFunction
        k: the diffusion coefficient
        f: the source term
        n: should be FacetNormal(mesh)
        penalty: the penalization of discontinuities and Dirichlet BCs
        dirichlet_bcs: a list of OcellarisDirichletBC objects
        neumann_bcs: a list of OcellarisNeumannBC objects
    
    Returns:
        a, L: the bilinear and linear forms
    """
    family = u.element().family()
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        eq = k*dot(nabla_grad(v), nabla_grad(u))*dx
        eq -= v*f*dx
        
        # Enforce Neumann BCs weakly
        for nbc in neumann_bcs:
            eq -= k*v*nbc.func()*nbc.ds()
    
    elif family == 'Discontinuous Lagrange':
        # Discontinous Galerkin implementation (Symmetric Interior Penalty method)
        eq = k*dot(nabla_grad(u), nabla_grad(v))*dx
        eq -= avg(k*dot(n, nabla_grad(u)))*jump(v)*dS
        eq -= avg(k*dot(n, nabla_grad(v)))*jump(u)*dS
        eq += avg(penalty)*jump(u)*jump(v)*dS
        eq -= f*v*dx
        
        # Enforce Dirichlet BCs weakly
        ds_penalty = dolfin.Constant(penalty*2)
        for dbc in dirichlet_bcs:
            eq -= k*dot(n, nabla_grad(u))*v*dbc.ds()
            eq -= k*dot(n, nabla_grad(v))*u*dbc.ds()
            eq += ds_penalty*u*v*dbc.ds()
            eq -= dbc.func()*(ds_penalty*v - k*dot(nabla_grad(v), n))*dbc.ds()
        
        # Enforce Neumann BCs weakly
        for nbc in neumann_bcs:
            eq -= nbc.func()*v*nbc.ds()
    
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)
    return a, L

def define_penalty(mesh, P, k_max, k_min):
    """
    Define the penalty parameter used in the Poisson equations
    
    Arguments:
        mesh: the mesh used in the simulation
        P: the polynomial degree of the unknown
        k_max: the maximum diffusion coefficient
        k_min: the minimum diffusion coefficient
    """
    ndim = mesh.geometry().dim()
    
    # Calculate geometrical factor used in the penalty
    geom_fac = 0
    for cell in dolfin.cells(mesh):
        vol = cell.volume()
        area = sum(cell.facet_area(i) for i in range(ndim + 1))
        gf = area/vol
        geom_fac = max(geom_fac, gf)
    geom_fac *= 1.0
    
    penalty = 300 * k_max**2/k_min * (P + 1)*(P + ndim)/ndim * geom_fac
    return penalty
