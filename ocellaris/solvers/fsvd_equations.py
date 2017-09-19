import dolfin
from dolfin import dot, grad, avg, jump, dx, dS, Constant
from . import BaseEquation
from ..solver_parts import define_penalty


class MomentumPredictionEquation(BaseEquation):
    def __init__(self, simulation, component, use_grad_p_form):
        """
        This class assembles the momentum equation for one velocity component, both CG and DG 
        """
        self.simulation = simulation
        self.component = component
        self.use_grad_p_form = use_grad_p_form
        self.define_momentum_equation()
        
    def calculate_penalties(self):
        """
        Calculate SIPG penalty
        """
        mpm = self.simulation.multi_phase_model
        mesh = self.simulation.data['mesh']
        
        mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
        P = self.simulation.data['Vu'].ufl_element().degree()
        penalty_dS = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0)
        penalty_ds = penalty_dS*2
        self.simulation.log.info('    DG SIP penalty viscosity:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
        D12 = Constant([0, 0])
        
        return Constant(penalty_dS), Constant(penalty_ds), D12
    
    def define_momentum_equation(self):
        """
        Setup the momentum equation for one velocity component
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mpm = self.simulation.multi_phase_model
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        ni = n[self.component]
        
        # Trial and test functions
        Vu = sim.data['Vu']
        u = dolfin.TrialFunction(Vu)
        v = dolfin.TestFunction(Vu)
        
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        u_conv = sim.data['u_conv']
        p = sim.data['p']
        
        # Fluid properties
        rho = mpm.get_density(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)
        
        # Values at previous time steps
        up = sim.data['up%d' % self.component]
        upp = sim.data['upp%d' % self.component]
        
        # Upwind and downwind velocities
        w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
        w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Penalties
        penalty_dS, penalty_ds, D12 = self.calculate_penalties()
        
        # Time derivative
        # ∂(ρu)/∂t
        eq = rho*(c1*u + c2*up + c3*upp)/dt*v*dx
        
        # Convection:
        # w⋅∇(ρu)    
        flux_nU = u*w_nU
        flux = jump(flux_nU)
        eq -= u*dot(grad(rho*v), u_conv)*dx
        eq += flux*jump(rho*v)*dS
        
        # Diffusion:
        # -∇⋅∇u
        eq += mu*dot(grad(u), grad(v))*dx
        
        # Symmetric Interior Penalty method for -∇⋅μ∇u
        eq -= avg(mu)*dot(n('+'), avg(grad(u)))*jump(v)*dS
        eq -= avg(mu)*dot(n('+'), avg(grad(v)))*jump(u)*dS
        
        # Symmetric Interior Penalty coercivity term
        eq += penalty_dS*jump(u)*jump(v)*dS
        
        # Pressure
        # ∇p
        if self.use_grad_p_form:
            eq += v*p.dx(self.component)*dx
            eq -= (avg(v) + dot(D12, jump(v, n)))*jump(p)*ni('+')*dS
        else:
            eq -= p*v.dx(self.component)*dx
            eq += (avg(p) - dot(D12, jump(p, n)))*jump(v)*ni('+')*dS
        
        # Body force (gravity)
        # ρ g
        eq -= rho*g[self.component]*v*dx
        
        # Other sources
        for f in sim.data['momentum_sources']:
            eq -= f[self.component]*v*dx
        
        # Dirichlet boundary
        dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % self.component, [])
        for dbc in dirichlet_bcs:
            u_bc = dbc.func()
            
            # Convection
            eq += rho*u*w_nU*v*dbc.ds()
            eq += rho*u_bc*w_nD*v*dbc.ds()
            
            # SIPG for -∇⋅μ∇u
            eq -= mu*dot(n, grad(u))*v*dbc.ds()
            eq -= mu*dot(n, grad(v))*u*dbc.ds()
            eq += mu*dot(n, grad(v))*u_bc*dbc.ds()
            
            # Weak Dirichlet
            eq += penalty_ds*u*v*dbc.ds()
            eq -= penalty_ds*u_bc*v*dbc.ds()
            
            # Pressure
            if not self.use_grad_p_form:
                eq += p*v*ni*dbc.ds()
        
        # Neumann boundary
        neumann_bcs = sim.data['neumann_bcs'].get('u%d' % self.component, [])
        for nbc in neumann_bcs:
            eq += mu*nbc.func()*v*nbc.ds()
            
            if not self.use_grad_p_form:
                eq += p*v*ni*nbc.ds()
        
        self.form_lhs, self.form_rhs = dolfin.system(eq)


class PressureCorrectionEquation(BaseEquation):
    def __init__(self, simulation):
        """
        This class assembles the pressure Poisson equation, both CG and DG 
        """
        self.simulation = simulation
        self.define_pressure_equation()
        
    def calculate_penalties(self):
        """
        Calculate SIPG penalty
        """
        mesh = self.simulation.data['mesh']
        P = self.simulation.data['Vp'].ufl_element().degree()
        k_min = k_max = 1.0
        penalty_dS = define_penalty(mesh, P, k_min, k_max, boost_factor=3, exponent=1.0)
        penalty_ds = penalty_dS*2
        self.simulation.log.info('    DG SIP penalty pressure:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
        return Constant(penalty_dS), Constant(penalty_ds)
    
    def define_pressure_equation(self):
        """
        Setup the pressure Poisson equation
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        Vp = sim.data['Vp']
        p_star = sim.data['p']
        u_star = sim.data['u']
        
        # Trial and test functions
        p = dolfin.TrialFunction(Vp)
        q = dolfin.TestFunction(Vp)
        
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties
        mpm = sim.multi_phase_model
        mu = mpm.get_laminar_dynamic_viscosity(0)
        rho_min = Constant(mpm.get_density_range()[0])
        
        # Weak form of the Poisson eq. with discontinuous elements
        # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^*
        K = 1.0
        a = K*dot(grad(p), grad(q))*dx
        L = K*dot(grad(p_star), grad(q))*dx

        # RHS, ∇⋅u^*
        u_flux = avg(u_star)
        L += c1*rho_min/dt*dot(u_star, grad(q))*dx
        L -= c1*rho_min/dt*dot(u_flux, n('+'))*jump(q)*dS
        
        # Symmetric Interior Penalty method for -∇⋅∇p
        a -= dot(n('+'), avg(K*grad(p)))*jump(q)*dS
        a -= dot(n('+'), avg(K*grad(q)))*jump(p)*dS
        
        # Symmetric Interior Penalty method for -∇⋅∇p^*
        L -= dot(n('+'), avg(K*grad(p_star)))*jump(q)*dS
        L -= dot(n('+'), avg(K*grad(q)))*jump(p_star)*dS
        
        # Weak continuity
        penalty_dS, penalty_ds = self.calculate_penalties()
        
        # Symmetric Interior Penalty coercivity term
        a += penalty_dS*jump(p)*jump(q)*dS
        #L += penalty_dS*jump(p_star)*jump(q)*dS
        
        # Collect Dirichlet and outlet boundary values
        dirichlet_vals_and_ds = []
        for dbc in sim.data['dirichlet_bcs'].get('p', []):
            dirichlet_vals_and_ds.append((dbc.func(), dbc.ds()))
        for obc in sim.data['outlet_bcs']:
            p_ = mu*dot(dot(grad(u_star), n), n)
            dirichlet_vals_and_ds.append((p_, obc.ds()))
        
        # Apply Dirichlet boundary conditions
        for p_bc, dds in dirichlet_vals_and_ds:
            # SIPG for -∇⋅∇p
            a -= dot(n, K*grad(p))*q*dds
            a -= dot(n, K*grad(q))*p*dds
            L -= dot(n, K*grad(q))*p_bc*dds
            
            # SIPG for -∇⋅∇p^*
            L -= dot(n, K*grad(p_star))*q*dds
            L -= dot(n, K*grad(q))*p_star*dds
            
            # Weak Dirichlet
            a += penalty_ds*p*q*dds
            L += penalty_ds*p_bc*q*dds
            
            # Weak Dirichlet for p^*
            #L += penalty_ds*p_star*q*dds
            #L -= penalty_ds*p_bc*q*dds
        
        # Neumann boundary conditions on p and p_star cancel
        #neumann_bcs = sim.data['neumann_bcs'].get('p', [])
        #for nbc in neumann_bcs:
        #    L += (nbc.func() - dot(n, grad(p_star)))*q*nbc.ds()
        
        # Use boundary conditions for the velocity for the
        # term from integration by parts of div(u_star)
        for d in range(sim.ndim):
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                L -= c1*rho_min/dt*u_bc*n[d]*q*dbc.ds()
            for nbc in neumann_bcs:
                L -= c1*rho_min/dt*u_star[d]*n[d]*q*nbc.ds()
        
        self.form_lhs = a
        self.form_rhs = L


class VelocityUpdateEquation(BaseEquation):
    def __init__(self, simulation, component):
        """
        Define the velocity update equation for velocity component d.
        """
        self.simulation = simulation
        self.component = component
        
        # Discontinuous or continuous elements
        Vu_family = simulation.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_update_equation()
    
    def define_update_equation(self):
        sim = self.simulation
        rho_min = Constant(self.simulation.multi_phase_model.get_density_range()[0])
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        
        Vu = sim.data['Vu']
        u_star = sim.data['u%d' % self.component]
        p_hat = sim.data['p_hat']
        u = dolfin.TrialFunction(Vu)
        v = dolfin.TestFunction(Vu)
        
        self.form_lhs = u*v*dx
        self.form_rhs = u_star*v*dx - dt/(c1*rho_min)*p_hat.dx(self.component)*v*dx


EQUATION_SUBTYPES = {
    'Default': (MomentumPredictionEquation, PressureCorrectionEquation, VelocityUpdateEquation),
}
