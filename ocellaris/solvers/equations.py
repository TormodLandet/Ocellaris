# encoding: utf8
import dolfin
from dolfin import dot, nabla_grad, avg, jump, dx, dS

def define_advection_problem(V, up, upp, u_conv, n, beta, time_coeffs, dt, dirichlet_bcs):
    """
    Define the advection problem
    
     d/dt(u) + u_conv ⋅ grad(u) = 0
     
    Returns the bilinear and linear forms
    """
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    
    # Upstream and downstream normal velocities
    flux_nU = u*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
    flux_nD = u*(dot(u_conv, n) - abs(dot(u_conv, n)))/2
    
    # Define the blended flux
    # The blending factor beta is not DG, so beta('+') == beta('-')
    b = beta('+')
    flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
    
    # Equation to solve
    c1, c2, c3 = time_coeffs 
    eq = (c1*u + c2*up + c3*upp)/dt*v*dx \
         - u*dot(u_conv, nabla_grad(v))*dx \
         + flux*jump(v)*dS
    
    # Enforce Dirichlet BCs weakly
    for dbc in dirichlet_bcs:
        eq += dot(u_conv, n)*v*(u - dbc.func())*dbc.ds()
    
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)
    
    return a, L

def define_poisson_problem(V, k, f, n, penalty, dirichlet_bcs, neumann_bcs):
    """
    Define the Poisson problem for u in f.space V
    
        - div(k*grad(u)) = f
    
    Note the minus in front of the first term!
    
    Returns the bilinear and linear forms
    """
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    
    # Interior
    a = dot(nabla_grad(v), nabla_grad(u))*dx
    
    # Interior facets
    a += avg(penalty)*jump(u)*jump(v)*dS \
         - dot(avg(nabla_grad(v)), n('+'))*jump(u)*dS \
         - dot(avg(nabla_grad(u)), n('+'))*jump(v)*dS

    # Source term
    L = -v*f/k*dx
    
    # Enforce Dirichlet BCs weakly
    for dbc in dirichlet_bcs:
        # These terms penalises jumps in (u - uD) on ds
        # where uD is the Dirichlet value on the boundary
        a += penalty*u*v*dbc.ds() \
             - dot(nabla_grad(v), n)*u*dbc.ds() \
             - dot(nabla_grad(u), n)*v*dbc.ds()
        L += penalty*v*dbc.func()*dbc.ds() \
             - dot(nabla_grad(v), n)*dbc.func()*dbc.ds()
    
    # Enforce Neumann BCs weakly
    for nbc in neumann_bcs:
        L += v*nbc.func()*nbc.ds()
    
    # FIXME: introduce k in the equation properly!!!!
    a = k*a
    L = k*L
    
    return a, L
