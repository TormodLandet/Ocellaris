# encoding: utf8
import dolfin
from dolfin import dot, nabla_grad, avg, jump, dx, ds, dS

def define_advection_problem(V, up, upp, u_conv, n, beta, dt):
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
    eq = (3*u - 4*up + upp)/(2*dt)*v*dx \
         - u*dot(u_conv, nabla_grad(v))*dx \
         + flux*jump(v)*dS \
         + u*dot(u_conv, n)*v*ds
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)
    
    return a, L

def define_poisson_problem(V, k, f, n, penalty, neumann_bcs):
    """
    Define the Poisson problem for u in f.space V
    
        div(k*grad(u)) = f
    
    Returns the bilinear and linear forms
    """
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    
    a = dot(nabla_grad(v), nabla_grad(u))*dx \
        - dot(avg(nabla_grad(v)), jump(u, n))*dS \
        - dot(jump(v, n), avg(nabla_grad(u)))*dS \
        + penalty*dot(jump(v, n), jump(u, n))*dS \
        - dot(nabla_grad(v), u*n)*ds \
        - dot(v*n, nabla_grad(u))*ds \
        + penalty*v*u*ds
    L = v*f*dx #- u0*dot(grad(v), n)*ds + (gamma/h)*u0*v*ds(1)
    
    # FIXME: introduce k in the equation properly!
    a = k*a
    
    # Add Neumann boundary conditions
    for nbc in neumann_bcs:
        L += v*nbc.value*nbc.ds 

    return a, L
