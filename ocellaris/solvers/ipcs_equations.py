# encoding: utf8
import dolfin
from dolfin import dot, nabla_grad, nabla_div, avg, jump, dx, dS
from . import BDF, CRANK_NICOLSON


# Default values, can be changed in the input file
PENALTY_BOOST = 2
STRESS_DIVERGENCE = True


class MomentumPredictionEquation(object):
    def __init__(self, simulation, d, beta, timestepping_method):
        """
        Define the momentum prediction equation for velocity component d.
        For DG "beta" is the convection blending factor (not used for CG)
        """
        self.simulation = simulation
        mesh = simulation.data['mesh']
        n = dolfin.FacetNormal(mesh)
        time_coeffs = simulation.data['time_coeffs']
        dt = simulation.data['dt']
        
        Vu = simulation.data['Vu']
        u_conv = simulation.data['u_conv']
        up = simulation.data['up%d' % d]
        upp = simulation.data['upp%d' % d]
        p = simulation.data['p']
        
        rho = simulation.data['rho']
        nu = simulation.data['nu']
        g = simulation.data['g']
        mu = rho*nu
        
        if up.element().family() == 'Discontinuous Lagrange':
            # Bounds on the diffusion coefficients for SIPG penalty calculation
            mpm = simulation.multi_phase_model
            mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
            k_u_max, k_u_min = mu_max, mu_min
            
            # Calculate SIPG penalties for the DG Poisson equation solver
            P = Vu.ufl_element().degree()
            boost_factor = simulation.input.get_value('solver/u/penalty_boost_factor', PENALTY_BOOST, 'float')
            penalty = define_penalty(mesh, P, k_u_max, k_u_min, boost_factor)
            simulation.log.info('\nDG SIPG penalties u%d:  %.4e' % (d, penalty))
            penalty = dolfin.Constant(penalty)
        else:
            penalty = None
        
        trial = dolfin.TrialFunction(Vu)
        test = dolfin.TestFunction(Vu)
        
        fp = -p.dx(d) + rho*g[d]
        dirichlet_bcs = simulation.data['dirichlet_bcs'].get('u%d' % d, [])
        neumann_bcs = simulation.data['neumann_bcs'].get('u%d' % d, [])
        
        if simulation.input.get_value('solver/use_stress_divergence_form', STRESS_DIVERGENCE, 'bool'):
            # Use the stress divergence form of the diffusion term
            #   Using u_conv here is unstable with BDF (oscillations)
            #   Using the previous time step values gives 2nd order
            #   convergence on the periodic Taylor-Green test case
            diff_u_expl = simulation.data['up'].dx(d)
        else:
            diff_u_expl = None
        
        if timestepping_method == BDF:
            thetas = dolfin.Constant([1.0, 0.0, 0.0])
            q = None
        elif timestepping_method == CRANK_NICOLSON:
            thetas = dolfin.Constant([0.5, 0.5, 0.0])
            q = up
        
        # Define:   ∂/∂t(ρ u) +  ∇⋅(ρ u ⊗ u_conv) = f_a = 0 + CN-terms
        a1, L1 = define_advection_problem(trial, test, up, upp, u_conv, rho, n, beta,
                                          time_coeffs, thetas, dt, dirichlet_bcs)
        
        # Define:   -∇⋅μ[(∇u) + (∇u_expl)^T] = f_p = - ∇p + ρg  + CN-terms 
        a2, L2 = define_poisson_problem(trial, test, mu, fp, n, thetas, q, penalty,
                                        dirichlet_bcs, neumann_bcs, diff_u_expl=diff_u_expl)
        
        self.form_lhs = a1+a2
        self.form_rhs = L1+L2
    
    def assemble_lhs(self):
        return dolfin.assemble(self.form_lhs)

    def assemble_rhs(self):
        return dolfin.assemble(self.form_rhs)


class PressureCorrectionEquation(object):
    def __init__(self, simulation):
        """
        Define the pressure correction equation
        """
        self.simulation = simulation
        mesh = simulation.data['mesh']
        n = dolfin.FacetNormal(mesh)
        time_coeffs = simulation.data['time_coeffs']
        dt = simulation.data['dt']
        rho = simulation.data['rho']
        
        Vu = simulation.data['Vu']
        Vp = simulation.data['Vp']
        u_star = simulation.data['u_star']
        p = simulation.data['p']
        
        if p.element().family() == 'Discontinuous Lagrange':
            # Bounds on the diffusion coefficients for SIPG penalty calculation
            mpm = simulation.multi_phase_model
            rho_max, rho_min = mpm.get_density_range()
            k_p_max, k_p_min = 1/rho_min, 1/rho_max
            
            # Calculate SIPG penalties for the DG Poisson equation solver
            Pu = Vu.ufl_element().degree()
            Pp = Vp.ufl_element().degree()
            P = max(Pu, Pp)
            boost_factor = simulation.input.get_value('solver/p/penalty_boost_factor', PENALTY_BOOST, 'float')
            penalty = define_penalty(mesh, P, k_p_max, k_p_min, boost_factor)
            simulation.log.info('\nDG SIPG penalties p:\n   %.4e' % penalty)
            penalty = dolfin.Constant(penalty)
        else:
            penalty = None
        
        trial = dolfin.TrialFunction(Vp)
        test = dolfin.TestFunction(Vp)
        f = -time_coeffs[0]/dt * nabla_div(u_star)
        k = 1/rho
        thetas = dolfin.Constant([1.0, -1.0])
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        neumann_bcs = self.simulation.data['neumann_bcs'].get('p', [])
        
        # Define:   -∇⋅[1/ρ (∇p_hat)] = - γ_1/Δt ∇⋅u_star    
        a, L = define_poisson_problem(trial, test, k, f, n, thetas, p,
                                      penalty, dirichlet_bcs, neumann_bcs)
        self.form_lhs = a
        self.form_rhs = L
    
    def assemble_lhs(self):
        return dolfin.assemble(self.form_lhs)

    def assemble_rhs(self):
        return dolfin.assemble(self.form_rhs)


class VelocityUpdateEquation(object):
    def __init__(self, simulation, d):
        """
        Define the velocity update equation for velocity component d.
        """
        self.simulation = simulation
        
        rho = simulation.data['rho']
        c0 = simulation.data['time_coeffs'][0]
        dt = simulation.data['dt']
        
        Vu = simulation.data['Vu']
        us = simulation.data['u_star%d' % d]
        p_hat = simulation.data['p_hat']
        
        u = dolfin.TrialFunction(Vu)
        v = dolfin.TestFunction(Vu)
        
        self.form_lhs = u*v*dx
        self.form_rhs = us*v*dx - dt/(c0*rho)*p_hat.dx(d)*v*dx
    
    def assemble_lhs(self):
        return dolfin.assemble(self.form_lhs)

    def assemble_rhs(self):
        return dolfin.assemble(self.form_rhs)


def define_advection_problem(u, v, up, upp, u_conv, r, n, beta, time_coeffs, thetas, dt, dirichlet_bcs):
    """
    Define the advection problem
    
        ∂/∂t(r u) +  ∇⋅(r U u_conv) = 0
    
    Where
    
        ∂u/∂t = γ1 u + γ2 up + γ3 upp
        γ1, γ2, γ3 = time_coeffs
        
        U = ϴ1 u + ϴ2 up + ϴ3 upp
        ϴ1, ϴ2, ϴ3 = thetas
        
        u = velocity component
        u_conv = velocity vector
        r = scalar (density)
    
    Returns the bilinear and linear forms
    """
    family = u.element().family()
    
    t1, t2, t3 = thetas
    U = t1*u + t2*up + t3*upp
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        c1, c2, c3 = time_coeffs
        eq = r*(c1*u + c2*up + c3*upp)/dt*v*dx + nabla_div(r*U*u_conv)*v*dx
    
    elif family == 'Discontinuous Lagrange':
        # Upstream and downstream normal velocities
        flux_nU = r*U*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
        flux_nD = r*U*(dot(u_conv, n) - abs(dot(u_conv, n)))/2
        
        # Define the blended flux
        # The blending factor beta is not DG, so beta('+') == beta('-')
        b = beta('+')
        flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
        
        # Equaqtion to solve
        c1, c2, c3 = time_coeffs 
        eq = r*(c1*u + c2*up + c3*upp)/dt*v*dx \
             - dot(r*U*u_conv, nabla_grad(v))*dx \
             + flux*jump(v)*dS
        
        # Enforce Dirichlet BCs weakly
        for dbc in dirichlet_bcs:
            eq += r*v*(u - dbc.func())*dbc.ds()
    
    a, L = dolfin.system(eq)    
    return a, L


def define_poisson_problem(u, v, k, f, n, thetas, q, penalty, dirichlet_bcs, neumann_bcs, diff_u_expl=None):
    """
    Define the Poisson problem for u:
    
        - ∇⋅(k∇U) = f
    
    or (if given vector diff_u_expl = ∂(u_expl)/∂x_i): 
    
        -∇⋅k[(∇U) + (∇u_expl)^T] = f
    
    where
    
        U = ϴ1 u + ϴ2 q
    
    Note the minus in front of the laplacian!
    
    Arguments:
        u: a TrialFunction
        v: a TestFunction
        k: the diffusion coefficient
        f: the source term
        n: should be FacetNormal(mesh)
        q: a known function (with the *same boundary conditions* as u)
        thetas: constant factors
        penalty: the penalization of discontinuities and Dirichlet BCs
        dirichlet_bcs: a list of OcellarisDirichletBC objects
        neumann_bcs: a list of OcellarisNeumannBC objects
        diff_u_expl: The velocity to use in the additional term in the
            stress divergence form of the viscous term in the momentum
            equations. Should be either u_expl.dx(i) or None
    
    Returns:
        a, L: the bilinear and linear forms
    """
    family = u.element().family()
    t1, t2 = thetas[0], thetas[1]
    
    if q is None:
        q = dolfin.Constant(0.0)
    
    U = t1*u + t2*q
    
    # Continuous Galerkin implementation of the Laplacian 
    eq = t1*k*dot(nabla_grad(U), nabla_grad(v))*dx - v*f*dx
    
    # Enforce Neumann BCs weakly. These are canceled by ∂q/∂n if ϴ1 + ϴ2 = 0
    for nbc in neumann_bcs:
        eq -= (t1 + t2)*k*nbc.func()*v*nbc.ds()
        
    # Extra terms in rhs due to special treatment of diffusive term in mom.eq
    if diff_u_expl is not None:
        eq += k*dot(diff_u_expl, nabla_grad(v))*dx
        for nbc in neumann_bcs:
            pass # TODO: include the boundary terms on the Neumann interfaces!
        
    if family == 'Discontinuous Lagrange':
        # Discontinous Galerkin implementation additional interior jump 
        # terms for the Symmetric Interior Penalty method
        eq -= avg(k*dot(n, nabla_grad(U)))*jump(v)*dS
        eq -= avg(k*dot(n, nabla_grad(v)))*jump(U)*dS
        eq += avg(penalty)*jump(u)*jump(v)*dS
        
        # Enforce Dirichlet BCs weakly
        ds_penalty = dolfin.Constant(penalty*2)
        for dbc in dirichlet_bcs:
            eq -= k*dot(n, nabla_grad(u))*v*dbc.ds()
            eq -= k*dot(n, nabla_grad(v))*u*dbc.ds()
            eq += ds_penalty*(u - dbc.func())*v*dbc.ds()
            eq += k*dot(nabla_grad(v), n)*dbc.func()*dbc.ds()
    
    a, L = dolfin.system(eq)
    return a, L


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
    for cell in dolfin.cells(mesh):
        vol = cell.volume()
        area = sum(cell.facet_area(i) for i in range(ndim + 1))
        gf = area/vol
        geom_fac = max(geom_fac, gf)
    
    penalty = boost_factor * k_max**2/k_min * (P + 1)*(P + ndim)/ndim * geom_fac**exponent
    return penalty
