# encoding: utf8
import dolfin
from dolfin import dot, div, grad, avg, jump, dx, dS, Constant
from . import BDF, CRANK_NICOLSON, UPWIND
from .dg_helpers import define_penalty


class MomentumPredictionEquation(object):
    def __init__(self, simulation, component, timestepping_method, flux_type,
                       use_stress_divergence_form, use_grad_p_form):
        """
        Define the momentum prediction equation for velocity component d.
        For DG "beta" is the convection blending factor (not used for CG)
        """
        self.simulation = simulation
        d = component
        self.timestepping_method = timestepping_method
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.flux_type = flux_type
        
        beta = Constant(0.0)
        mesh = simulation.data['mesh']
        n = dolfin.FacetNormal(mesh)
        time_coeffs = simulation.data['time_coeffs']
        dt = simulation.data['dt']
        
        Vu = simulation.data['Vu']
        u_conv = simulation.data['u_conv']
        up = simulation.data['up%d' % d]
        upp = simulation.data['upp%d' % d]
        p = simulation.data['p']
        
        rho = simulation.data['rho_star']
        nu = simulation.data['nu_star']
        g = simulation.data['g']
        mu = rho*nu
        
        if up.ufl_element().family() == 'Discontinuous Lagrange':
            # Bounds on the diffusion coefficients for SIPG penalty calculation
            mpm = simulation.multi_phase_model
            mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
            k_u_max, k_u_min = mu_max, mu_min
            
            # Calculate SIPG penalties for the DG Poisson equation solver
            P = Vu.ufl_element().degree()
            penalty = define_penalty(mesh, P, k_u_max, k_u_min)
            simulation.log.info('\nDG SIPG penalties u%d:  %.4e' % (d, penalty))
            penalty = dolfin.Constant(penalty)
        else:
            penalty = None
        
        trial = dolfin.TrialFunction(Vu)
        test = dolfin.TestFunction(Vu)
        
        fp = -p.dx(d) + rho*g[d]
        dirichlet_bcs = simulation.data['dirichlet_bcs'].get('u%d' % d, [])
        neumann_bcs = simulation.data['neumann_bcs'].get('u%d' % d, [])
        
        if self.use_stress_divergence_form:
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
    def __init__(self, simulation, use_lagrange_multiplicator):
        """
        Define the pressure correction equation
        """
        self.simulation = simulation
        mesh = simulation.data['mesh']
        n = dolfin.FacetNormal(mesh)
        time_coeffs = simulation.data['time_coeffs']
        dt = simulation.data['dt']
        rho = simulation.data['rho_star']
        
        Vu = simulation.data['Vu']
        Vp = simulation.data['Vp']
        u_star = simulation.data['u_star']
        p = simulation.data['p']
        
        if p.ufl_element().family() == 'Discontinuous Lagrange':
            # Bounds on the diffusion coefficients for SIPG penalty calculation
            #mpm = simulation.multi_phase_model
            #rho_max, rho_min = mpm.get_density_range()
            #k_p_max, k_p_min = 1/rho_min, 1/rho_max
            
            # Calculate SIPG penalties for the DG Poisson equation solver
            Pu = Vu.ufl_element().degree()
            Pp = Vp.ufl_element().degree()
            P = max(Pu, Pp)
            penalty = define_penalty(mesh, P, 1, 1)
            simulation.log.info('\nDG SIPG penalties p:\n   %.4e' % penalty)
            penalty = dolfin.Constant(penalty)
        else:
            penalty = None
        
        trial = dolfin.TrialFunction(Vp)
        test = dolfin.TestFunction(Vp)
        f = -rho*time_coeffs[0]/dt * div(u_star)
        k = Constant(1)
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
        
        rho = simulation.data['rho_star']
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
        
        u = velocity componentcomponent
        u_conv = velocity vector
        r = scalar (density)
    
    Returns the bilinear and linear forms
    """
    family = u.ufl_element().family()
    
    t1, t2, t3 = thetas
    U = t1*u + t2*up + t3*upp
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        c1, c2, c3 = time_coeffs
        eq = r*(c1*u + c2*up + c3*upp)/dt*v*dx + r*div(U*u_conv)*v*dx
    
    elif family == 'Discontinuous Lagrange':
        # Upstream and downstream normal velocities
        flux_nU = r*U*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
        flux_nD = r*U*(dot(u_conv, n) - abs(dot(u_conv, n)))/2
        
        # Define the blended flux
        # The blending factor beta is not DG, so beta('+') == beta('-')
        b = beta('+')
        flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
        
        # Equation to solve
        c1, c2, c3 = time_coeffs 
        eq = r*(c1*u + c2*up + c3*upp)/dt*v*dx \
             - dot(r*U*u_conv, grad(v))*dx \
             + flux*jump(v)*dS
        
        # Enforce Dirichlet BCs weakly
        for dbc in dirichlet_bcs:
            eq += dot(u_conv, n)*r*v*(u - dbc.func())*dbc.ds()
    
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
        
        
        a2, L2 = define_poisson_problem(trial, test, mu, fp, n, thetas, q, penalty,
                                        dirichlet_bcs, neumann_bcs, diff_u_expl=diff_u_expl)
                                        
        a, L = define_poisson_problem(trial, test, k, f, n, thetas, p,
                                      penalty, dirichlet_bcs, neumann_bcs)
    """
    family = u.ufl_element().family()
    t1, t2 = thetas[0], thetas[1]
    
    if q is None:
        q = dolfin.Constant(0.0)
    
    U = t1*u + t2*q
    
    # Continuous Galerkin implementation of the Laplacian 
    eq = t1*k*dot(grad(U), grad(v))*dx - v*f*dx
    
    # Enforce Neumann BCs weakly. These are canceled by ∂q/∂n if ϴ1 + ϴ2 = 0
    for nbc in neumann_bcs:
        eq -= (t1 + t2)*k*nbc.func()*v*nbc.ds()
        
    # Extra terms in rhs due to special treatment of diffusive term in mom.eq
    if diff_u_expl is not None:
        eq += k*dot(diff_u_expl, grad(v))*dx
        for nbc in neumann_bcs:
            pass # TODO: include the boundary terms on the Neumann interfaces!
        
    if family == 'Discontinuous Lagrange':
        # Discontinous Galerkin implementation additional interior jump 
        # terms for the Symmetric Interior Penalty method
        eq -= avg(k*dot(n, grad(U)))*jump(v)*dS
        eq -= avg(k*dot(n, grad(v)))*jump(U)*dS
        eq += avg(penalty)*jump(u)*jump(v)*dS
        
        # Enforce Dirichlet BCs weakly
        ds_penalty = dolfin.Constant(penalty*2)
        for dbc in dirichlet_bcs:
            eq -= k*dot(n, grad(u))*v*dbc.ds()
            eq -= k*dot(n, grad(v))*u*dbc.ds()
            eq += ds_penalty*(u - dbc.func())*v*dbc.ds()
            eq += k*dot(grad(v), n)*dbc.func()*dbc.ds()
    
    a, L = dolfin.system(eq)
    return a, L


# ----------------------------------------------------------------------------------------------------------------------
class MomentumPredictionEquationNew(object):
    def __init__(self, simulation, component, timestepping_method, flux_type,
                 use_stress_divergence_form, use_grad_p_form):
        """
        This class assembles the momentum equation for one velocity component, both CG and DG 
        """
        self.simulation = simulation
        self.component = component
        self.timestepping_method = timestepping_method
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.flux_type = flux_type
        
        # Discontinuous or continuous elements
        Vu_family = simulation.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_momentum_equation()
        
    def calculate_penalties(self, nu):
        """
        Calculate SIPG penalty
        """
        mpm = self.simulation.multi_phase_model
        mesh = self.simulation.data['mesh']
        
        mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
        P = self.simulation.data['Vu'].ufl_element().degree()
        penalty_dS = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0)
        penalty_ds = penalty_dS*2
        self.simulation.log.info('DG SIP penalty viscosity:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
        D12 = Constant([1, 1])
        
        return Constant(penalty_dS), Constant(penalty_ds), D12
    
    def define_momentum_equation(self):
        """
        Setup the momentum equation for one velocity component
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Trial and test functions
        Vu = sim.data['Vu']
        u = dolfin.TrialFunction(Vu)
        v = dolfin.TestFunction(Vu)
        
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        u_conv = sim.data['u_conv']
        p = sim.data['p']
        
        # Fluid properties at t^{n+1}*
        rhos = sim.data['rho_star']
        nus = sim.data['nu_star']
        mus = rhos*nus
        
        # Include (∇u)^T term?
        assert not self.use_stress_divergence_form
        
        if self.vel_is_discontinuous:
            penalty_dS, penalty_ds, D12 = self.calculate_penalties(nus)
            
            # Upwind and downwind velocitues
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Values at previous time steps
        up = sim.data['up%d' % self.component]
        upp = sim.data['upp%d' % self.component]
            
        if not self.vel_is_discontinuous:
            # Weak form of the Navier-Stokes eq. with continuous elements
            
            # Time derivative
            # ∂u/∂t
            a = rhos*c1*u/dt*v*dx
            L = -rhos*(c2*up + c3*upp)/dt*v*dx
            
            # Convection
            # ρ∇⋅(u ⊗ u_conv)
            a += rhos*div(u*u_conv)*v*dx
            
            # Diffusion
            # -∇⋅μ∇u
            a += mus*dot(grad(u), grad(v))*dx
            
            # Pressure
            # ∇p
            if self.use_grad_p_form:
                L -= p.dx(self.component)*v*dx
            else:
                L += v.dx(self.component)*p*dx
            
            # Body force (gravity)
            # ρ g
            L += rhos*g[self.component]*v*dx
            
            # Neumann boundary
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % self.component, [])
            for nbc in neumann_bcs:
                L += mus*nbc.func()*v*nbc.ds()
                
                if not self.use_grad_p_form:
                    L -= p*v*n[self.component]*nbc.ds()
        
        else:
            # Weak form of the Navier-Stokes eq. with discontinuous elements
            assert self.flux_type == UPWIND
            
            # Time derivative
            # ∂u/∂t
            a = rhos*c1*u/dt*v*dx
            L = -rhos*(c2*up + c3*upp)/dt*v*dx
            
            # Convection:
            # -w⋅∇u    
            flux_nU = rhos*u*w_nU
            flux = jump(flux_nU)
            a -= rhos*u*div(v*u_conv)*dx
            a += flux*jump(v)*dS
            
            # Diffusion:
            # -∇⋅∇u
            a += mus*dot(grad(u), grad(v))*dx
            
            # Symmetric Interior Penalty method for -∇⋅μ∇u
            a -= avg(mus)*dot(n('+'), avg(grad(u)))*jump(v)*dS
            a -= avg(mus)*dot(n('+'), avg(grad(v)))*jump(u)*dS
            
            # Symmetric Interior Penalty coercivity term
            a += penalty_dS*jump(u)*jump(v)*dS
            
            # Pressure
            # ∇p
            if self.use_grad_p_form:
                L -= v*p.dx(self.component)*dx
                L += (avg(v) + D12*jump(v, n))*jump(p)*n[self.component]('+')*dS
            else:
                L += p*v.dx(self.component)*dx
                L -= (avg(p) - dot(D12, jump(p, n)))*jump(v)*n[self.component]('+')*dS
            
            # Body force (gravity)
            # ρ g
            L += rhos*g[self.component]*v*dx
            
            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % self.component, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                
                # Convection
                L -= rhos*w_nD*u_bc*v*dbc.ds()
                
                # SIPG for -∇⋅μ∇u
                a -= mus*dot(n, grad(u))*v*dbc.ds()
                a -= mus*dot(n, grad(v))*u*dbc.ds()
                L -= mus*dot(n, grad(v))*u_bc*dbc.ds()
                
                # Weak Dirichlet
                a += penalty_ds*u*v*dbc.ds()
                L += penalty_ds*u_bc*v*dbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    L -= p*v*n[self.component]*dbc.ds()
        
        self.form_lhs = a
        self.form_rhs = L
        self.tensor_lhs = None
        self.tensor_rhs = None
    
    def assemble_lhs(self):
        if self.tensor_lhs is None:
            self.tensor_lhs = dolfin.assemble(self.form_lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            self.tensor_rhs = dolfin.assemble(self.form_rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs


class PressureCorrectionEquationNew(object):
    def __init__(self, simulation, use_lagrange_multiplicator):
        """
        This class assembles the pressure Poisson equation, both CG and DG 
        """
        self.simulation = simulation
        self.use_lagrange_multiplicator = use_lagrange_multiplicator
        
        # Discontinuous or continuous elements
        Vp_family = simulation.data['Vp'].ufl_element().family()
        self.pressure_is_discontinuous = (Vp_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_pressure_equation()
        
    def calculate_penalties(self):
        """
        Calculate SIPG penalty
        """
        mesh = self.simulation.data['mesh']
        P = self.simulation.data['Vp'].ufl_element().degree()
        k_min = k_max = 1
        penalty_dS = define_penalty(mesh, P, k_min, k_max, boost_factor=3, exponent=1.0)
        penalty_ds = penalty_dS*2
        self.simulation.log.info('DG SIP penalty pressure:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
                
        return Constant(penalty_dS), Constant(penalty_ds)
    
    def define_pressure_equation(self):
        """
        Setup the pressure Poisson equation
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        Vp = sim.data['Vp']
        p_star = sim.data['p']
        u_star = sim.data['u_star']
        
        # Trial and test functions
        p = dolfin.TrialFunction(Vp)
        q = dolfin.TestFunction(Vp) 
        
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties at t^{n+1}*
        rhos = sim.data['rho_star']
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        assert not self.use_lagrange_multiplicator, 'NOT IMPLEMENTED YET'
        
        if not self.pressure_is_discontinuous:
            # Weak form of the Poisson eq. with continuous elements
            # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^* 
            a = dot(grad(p), grad(q))*dx
            L = dot(grad(p_star), grad(q))*dx
            L -= c1/dt*rhos*div(u_star)*q*dx
            
            # Neumann boundary conditions on p and p_star cancel
        
        else:
            # Weak form of the Poisson eq. with discontinuous elements
            # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^*
            a = dot(grad(p), grad(q))*dx
            L = dot(grad(p_star), grad(q))*dx
            L -= c1/dt*rhos*div(u_star)*q*dx
            
            # Symmetric Interior Penalty method for -∇⋅∇p
            a -= dot(n('+'), avg(grad(p)))*jump(q)*dS
            a -= dot(n('+'), avg(grad(q)))*jump(p)*dS
            
            # Symmetric Interior Penalty method for -∇⋅∇p^*
            L -= dot(n('+'), avg(grad(p_star)))*jump(q)*dS
            L -= dot(n('+'), avg(grad(q)))*jump(p_star)*dS
            
            # Weak continuity
            penalty_dS, penalty_ds = self.calculate_penalties()
            
            # Symmetric Interior Penalty coercivity term
            a += penalty_dS*jump(p)*jump(q)*dS
            
            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('p', [])
            for dbc in dirichlet_bcs:
                p_bc = dbc.func()
                
                # SIPG for -∇⋅∇p
                a -= dot(n, grad(p))*q*dbc.ds()
                a -= dot(n, grad(q))*p*dbc.ds()
                L -= dot(n, grad(q))*p_bc*dbc.ds()
                
                L -= dot(n, grad(p_star))*q*dbc.ds()
                L -= dot(n, grad(q))*p_star*dbc.ds()
                
                # Weak Dirichlet
                a += penalty_ds*p*q*dbc.ds()
                L += penalty_ds*p_bc*q*dbc.ds()
            
            # Neumann boundary conditions on p and p_star cancel
        
        self.form_lhs = a
        self.form_rhs = L
        self.tensor_lhs = None
        self.tensor_rhs = None
    
    def assemble_lhs(self):
        if self.tensor_lhs is None:
            self.tensor_lhs = dolfin.assemble(self.form_lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            self.tensor_rhs = dolfin.assemble(self.form_rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs

EQUATION_SUBTYPES = {
    'Old': (MomentumPredictionEquation, PressureCorrectionEquation, VelocityUpdateEquation),
    'New': (MomentumPredictionEquationNew, PressureCorrectionEquationNew, VelocityUpdateEquation)
}