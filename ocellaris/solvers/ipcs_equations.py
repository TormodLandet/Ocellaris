# encoding: utf8
import dolfin
from dolfin import dot, nabla_grad, nabla_div, avg, jump, dx, dS
from . import BDF, CRANK_NICOLSON


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
        
        if up.element().family() == 'Discontinuous Lagrange':
            # Bounds on the diffusion coefficients for SIPG penalty calculation
            mpm = simulation.multi_phase_model
            nu_max, nu_min = mpm.get_laminar_kinematic_viscosity_range()
            k_u_max, k_u_min = nu_max, nu_min
            
            # Calculate SIPG penalties for the DG Poisson equation solver
            P = Vu.ufl_element().degree()
            penalty = define_penalty(mesh, P, k_u_max, k_u_min)
            simulation.log.info('\nDG SIPG penalties u%d:  %.4e' % (d, penalty))
            penalty = dolfin.Constant(penalty)
        else:
            penalty = None
        
        trial = dolfin.TrialFunction(Vu)
        test = dolfin.TestFunction(Vu)
        
        fp = -1/rho*p.dx(d) + g[d]
        dirichlet_bcs = simulation.data['dirichlet_bcs'].get('u%d' % d, [])
        neumann_bcs = simulation.data['neumann_bcs'].get('u%d' % d, [])
        
        # Use the stress divergence form of the diffusion term
        if simulation.input.get_value('solver/use_stress_divergence_form', True, 'bool'):
            diff_u_expl = u_conv.dx(d)
        else:
            diff_u_expl = None
        
        if timestepping_method == BDF:
            uc = u_conv
            fa = dolfin.Constant(0)
        elif timestepping_method == CRANK_NICOLSON:
            theta = dolfin.Constant(0.5)
            uc = theta*u_conv
            fa = -dot(uc, nabla_grad(up[d]))
            fp = fp/theta + nabla_div(nu*nabla_grad(up))
        
        a1, L1 = define_advection_problem(trial, test, up, upp, uc, fa,
                                          n, beta, time_coeffs, dt, dirichlet_bcs)
        a2, L2 = define_poisson_problem(trial, test, nu, fp, n, penalty,
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
            penalty = define_penalty(mesh, P, k_p_max, k_p_min)
            simulation.log.info('\nDG SIPG penalties p:\n   p: %.4e' % penalty)
            penalty = dolfin.Constant(penalty)
        else:
            penalty = None
        
        trial = dolfin.TrialFunction(Vp)
        test = dolfin.TestFunction(Vp)
        f = -time_coeffs[0]/dt * nabla_div(u_star)
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        neumann_bcs = self.simulation.data['neumann_bcs'].get('p', [])
        a, L = define_poisson_problem(trial, test, 1/rho, f, n, penalty,
                                      dirichlet_bcs, neumann_bcs, q=p)
        self.form_lhs = a
        self.form_rhs = L
    
    def assemble_lhs(self):
        return dolfin.assemble(self.form_lhs)

    def assemble_rhs(self):
        return dolfin.assemble(self.form_rhs)


def define_advection_problem(u, v, up, upp, u_conv, f, n, beta, time_coeffs, dt, dirichlet_bcs):
    """
    Define the advection problem
    
     d/dt(u) + u_conv ⋅ grad(u) = f
     
    Returns the bilinear and linear forms
    """
    family = u.element().family()
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        c1, c2, c3 = time_coeffs 
        eq = (c1*u + c2*up + c3*upp)/dt*v*dx + dot(u_conv, nabla_grad(u))*v*dx - f*v*dx
    
    elif family == 'Discontinuous Lagrange':
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
             + flux*jump(v)*dS  - f*v*dx
        
        # Enforce Dirichlet BCs weakly
        for dbc in dirichlet_bcs:
            eq += dot(u_conv, n)*v*(u - dbc.func())*dbc.ds()
        
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)    
    return a, L


def define_poisson_problem(u, v, k, f, n, penalty, dirichlet_bcs, neumann_bcs, q=None, diff_u_expl=None):
    """
    Define the Poisson problem for u
    
        - div(k*grad(u - q)) = f
    
    Note the minus in front of the first term! The known value q does
    not have to be given. If it is it will be included in the linear
    form that is returned.
    
    Arguments:
        u: a TrialFunction
        v: a TestFunction
        k: the diffusion coefficient
        f: the source term
        n: should be FacetNormal(mesh)
        penalty: the penalization of discontinuities and Dirichlet BCs
        dirichlet_bcs: a list of OcellarisDirichletBC objects
        neumann_bcs: a list of OcellarisNeumannBC objects
        q: a known function
        diff_u_expl: The velocity to use in the additional term in the
            stress divergence form of the viscous term in the momentum
            equations. Should be either u_expl.dx(i) or None
    
    Returns:
        a, L: the bilinear and linear forms
    """
    family = u.element().family()
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        eq = k*dot(nabla_grad(v), nabla_grad(u))*dx
        eq -= v*f*dx
    
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
        eq -= nbc.func()*v*nbc.ds() # TODO: is this missing a "k"?
    
    # Extra terms in rhs due to q != 0       
    if q is not None:
        eq -= k*dot(nabla_grad(q), nabla_grad(v))*dx
        for nbc in neumann_bcs:
            eq += dot(nabla_grad(q), n)*v*nbc.ds()
    
    # Extra terms in rhs due to special treatment of diffusive term in mom.eq
    if diff_u_expl is not None:
        eq += k*dot(nabla_grad(v), diff_u_expl)*dx
        for nbc in neumann_bcs:
            pass # TODO: include the boundary terms on the Neumann interfaces!
    
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
