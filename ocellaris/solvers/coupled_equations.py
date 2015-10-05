# encoding: utf8
import dolfin
from dolfin import dx, div, grad, dot, jump, avg, ds, dS, Constant
from ufl.classes import Zero
from . import BDF, UPWIND
from .dg_helpers import define_penalty

# Default values, can be changed in the input file
PENALTY_BOOST = 1
STRESS_DIVERGENCE = True


class CoupledEquationsSlow(object):
    def __init__(self, simulation, timestepping_method):
        """
        Define the coupled Navier-Stokes equations

        This implementation assembles the full LHS and RHS each time they are needed
        """
        self.simulation = simulation
        self.timestepping_method = timestepping_method

        # Read necessary input
        inp = simulation.input

        # Get flux type for the velocity
        self.flux_type = inp.get_value('convection/u/flux_type', UPWIND, 'string')
        
        # Stress divergence form
        self.use_stress_divergence_form = inp.get_value('solver/use_stress_divergence_form', STRESS_DIVERGENCE, 'bool')

        # Discontinuous or continuous elements
        Vu_family = simulation.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_coupled_equation()
    
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        """
        sim = self.simulation
        Vcoupled = sim.data['Vcoupled']
        u_conv = sim.data['u_conv']
        
        # Unpack the coupled trial and test functions
        uc = dolfin.TrialFunction(Vcoupled)
        vc = dolfin.TestFunction(Vcoupled)
        ulist = []; vlist = []
        ndim = self.simulation.ndim
        for d in range(ndim):
            ulist.append(uc[d])
            vlist.append(vc[d])
        
        u = dolfin.as_vector(ulist)
        v = dolfin.as_vector(vlist)
        p = uc[ndim]
        q = vc[ndim]
        lm_trial = uc[ndim+1]
        lm_test = vc[ndim+1]
        
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties at t^{n}, t^{n-1} and t^{n+1}*
        rhop = sim.data['rho']
        rhopp = sim.data['rho_old']
        rhos = sim.data['rho_star']
        nus = sim.data['nu_star']
        mus = rhos*nus
        
        # Include (∇u)^T term?
        if self.use_stress_divergence_form:
            sd = Constant(1.0)
        else:
            sd = Constant(0.0)
            
        if self.vel_is_discontinuous:
            # Calculate SIPG penalty
            mpm = sim.multi_phase_model
            mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
            P = sim.data['Vu'].ufl_element().degree()
            penalty1 = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0)
            penalty2 = penalty1*2
            penalty_dS = Constant(penalty1)
            penalty_ds = Constant(penalty2)
            sim.log.info('DG SIP penalty:  dS %.1f  ds %.1f' % (penalty1, penalty2))
            
            # Upwind and downwind velocitues
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        eq = (p*lm_test + q*lm_trial)*dx
        
        # Momentum equations
        for d in range(sim.ndim):
            up = sim.data['up%d' % d]
            upp = sim.data['upp%d' % d]
            
            if not self.vel_is_discontinuous:
                # Weak form of the Navier-Stokes eq. with continuous elements
                
                # Divergence free criterion
                # ∇⋅u = 0
                eq += u[d].dx(d)*q*dx
                
                # Time derivative
                # ∂u/∂t
                eq += (rhos*c1*u[d] + rhop*c2*up + rhopp*c3*upp)/dt*v[d]*dx
                
                # Convection
                # ∇⋅(ρ u ⊗ u_conv)
                eq += div(rhos*u[d]*u_conv)*v[d]*dx
                
                # Diffusion
                # -∇⋅μ[(∇u) + (∇u)^T]
                eq += mus*dot(grad(u[d]), grad(v[d]))*dx
                eq += sd*mus*dot(u.dx(d), grad(v[d]))*dx
                
                # Pressure
                # ∇p
                eq -= v[d].dx(d)*p*dx
                
                # Body force (gravity)
                # ρ g
                eq -= rhos*g[d]*v[d]*dx
                
                # Neumann boundary conditions
                neumann_bcs_pressure = sim.data['neumann_bcs'].get('p', [])
                for nbc in neumann_bcs_pressure:
                    eq += p*v[d]*n[d]*nbc.ds()
                
            else:
                # Weak form of the Navier-Stokes eq. with discontinuous elements
                assert self.flux_type == UPWIND
                
                # Divergence free criterion
                # ∇⋅u = 0
                eq -= u[d]*q.dx(d)*dx
                eq += avg(u[d])*jump(q)*n[d]('+')*dS
                
                # Time derivative
                # ∂u/∂t
                eq += (rhos*c1*u[d] + rhop*c2*up + rhopp*c3*upp)/dt*v[d]*dx
                
                # Convection:
                # -w⋅∇u    
                flux_nU = rhos*u[d]*w_nU
                flux = jump(flux_nU)
                eq -= rhos*u[d]*div(v[d]*u_conv)*dx
                eq += flux*jump(v[d])*dS
                eq += rhos*flux_nU*v[d]*ds
                
                # Diffusion:
                # -∇⋅∇u
                assert not self.use_stress_divergence_form
                eq += mus*dot(grad(u[d]), grad(v[d]))*dx
                
                # Symmetric Interior Penalty method for -∇⋅μ∇u
                eq -= avg(mus)*dot(n('+'), avg(grad(u[d])))*jump(v[d])*dS
                eq -= avg(mus)*dot(n('+'), avg(grad(v[d])))*jump(u[d])*dS
                
                # Symmetric Interior Penalty coercivity term
                eq += penalty_dS*jump(u[d])*jump(v[d])*dS
                
                # Pressure
                # ∇p
                eq -= v[d].dx(d)*p*dx
                eq += avg(p)*jump(v[d])*n[d]('+')*dS
                
                # Body force (gravity)
                # ρ g
                eq -= rhos*g[d]*v[d]*dx
                
                # Dirichlet boundary
                dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
                for dbc in dirichlet_bcs:
                    # Divergence free criterion
                    eq += q*u[d]*n[d]*dbc.ds()
                    
                    # Convection
                    eq += rhos*w_nD*dbc.func()*v[d]*dbc.ds()
                    
                    # SIPG for -∇⋅μ∇u
                    eq -= mus*dot(n, grad(u[d]))*v[d]*dbc.ds()
                    eq -= mus*dot(n, grad(v[d]))*u[d]*dbc.ds()
                    eq += mus*dot(n, grad(v[d]))*dbc.func()*dbc.ds()
                    
                    # Weak Dirichlet
                    eq += penalty_ds*(u[d] - dbc.func())*v[d]*dbc.ds()
                    
                    # Pressure
                    eq += p*v[d]*n[d]*dbc.ds()
                
        a, L = dolfin.system(eq)
        self.form_lhs = a
        self.form_rhs = L
    
    def assemble_lhs(self):
        return dolfin.assemble(self.form_lhs)

    def assemble_rhs(self):
        return dolfin.assemble(self.form_rhs)


class CoupledEquationsFast(object):
    def __init__(self, simulation, timestepping_method):
        """
        Define the coupled Navier-Stokes equations
        
        This implementation tries to avoid re-assembly of matrices as much as possible
        """
        self.simulation = simulation
        self.timestepping_method = timestepping_method

        # Read necessary input
        inp = simulation.input

        # Get flux type for the velocity
        self.flux_type = inp.get_value('convection/u/flux_type', UPWIND, 'string')
        
        # Stress divergence form
        self.use_stress_divergence_form = inp.get_value('solver/use_stress_divergence_form', STRESS_DIVERGENCE, 'bool')

        # Discontinuous or continuous elements
        Vu_family = simulation.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_coupled_equation()
    
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        """
        assert self.vel_is_discontinuous and not self.use_stress_divergence_form
        assert self.flux_type == UPWIND
        
        sim = self.simulation
        Vcoupled = sim.data['Vcoupled']
        u_conv = sim.data['u_conv']
        
        # Unpack the coupled trial and test functions
        uc = dolfin.TrialFunction(Vcoupled)
        vc = dolfin.TestFunction(Vcoupled)
        ulist = []; vlist = []
        ndim = self.simulation.ndim
        for d in range(ndim):
            ulist.append(uc[d])
            vlist.append(vc[d])
        
        u = dolfin.as_vector(ulist)
        v = dolfin.as_vector(vlist)
        p = uc[ndim]
        q = vc[ndim]
        lm_trial = uc[ndim+1]
        lm_test = vc[ndim+1]
        
        _, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties at t^{n}, t^{n-1} and t^{n+1}*
        rhop = sim.data['rho']
        rhopp = sim.data['rho_old']
        rhos = sim.data['rho_star']
        nus = sim.data['nu_star']
        mus = rhos*nus
            
        if self.vel_is_discontinuous:
            # Calculate SIPG penalty
            mpm = sim.multi_phase_model
            mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
            P = sim.data['Vu'].ufl_element().degree()
            penalty1 = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0)
            penalty2 = penalty1*2
            penalty_dS = Constant(penalty1)
            penalty_ds = Constant(penalty2)
            sim.log.info('DG SIP penalty:  dS %.1f  ds %.1f' % (penalty1, penalty2))
            
            # Upwind and downwind velocitues
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        a_lm = (p*lm_test + q*lm_trial)*dx

        # Forms to be assembled for each dimension
        empty = Zero()*dx(domain=mesh)     
        a_divergence = empty
        a_convection = empty
        a_diffusion = empty
        a_jump = empty
        a_pressure = empty
        a_dirichlet = empty
        L = empty
        
        # Momentum equations
        for d in range(sim.ndim):
            up = sim.data['up%d' % d]
            upp = sim.data['upp%d' % d]
            
            # Divergence free criterion
            # ∇⋅u = 0
            a_divergence -= u[d]*q.dx(d)*dx
            a_divergence += avg(u[d])*jump(q)*n[d]('+')*dS
            
            # Time derivative
            # ∂u/∂t
            L -= (rhop*c2*up + rhopp*c3*upp)/dt*v[d]*dx
            
            # Convection:
            # -w⋅∇u    
            flux_nU = u[d]*w_nU
            flux = jump(flux_nU)
            a_convection -= rhos*u[d]*div(v[d]*u_conv)*dx
            a_convection += rhos*flux*jump(v[d])*dS
            a_convection += rhos*flux_nU*v[d]*ds
            
            # Diffusion:
            # -∇⋅∇u
            a_diffusion += dot(grad(u[d]), grad(v[d]))*dx
            
            # Symmetric Interior Penalty method for -∇⋅μ∇u
            a_diffusion -= dot(n('+'), avg(grad(u[d])))*jump(v[d])*dS
            a_diffusion -= dot(n('+'), avg(grad(v[d])))*jump(u[d])*dS
            
            # Symmetric Interior Penalty coercivity term
            a_jump += penalty_dS*jump(u[d])*jump(v[d])*dS
            
            # Pressure
            # ∇p
            a_pressure -= v[d].dx(d)*p*dx
            a_pressure += avg(p)*jump(v[d])*n[d]('+')*dS
            
            # Body force (gravity)
            # ρ g
            L += rhos*g[d]*v[d]*dx
            
            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                # Divergence free criterion
                a_divergence += q*u[d]*n[d]*dbc.ds()
                
                # Convection
                L -= rhos*w_nD*dbc.func()*v[d]*dbc.ds()
                
                # SIPG for -∇⋅μ∇u
                a_diffusion -= dot(n, grad(u[d]))*v[d]*dbc.ds()
                a_diffusion -= dot(n, grad(v[d]))*u[d]*dbc.ds()
                L -= mus*dot(n, grad(v[d]))*dbc.func()*dbc.ds()
                
                # Weak Dirichlet
                a_dirichlet += penalty_ds*u[d]*v[d]*dbc.ds()
                L += penalty_ds*dbc.func()*v[d]*dbc.ds()
                
                # Pressure
                a_pressure += p*v[d]*n[d]*dbc.ds()

        # Assemble matrices
        self.matrix_Q = dolfin.assemble(a_lm + a_divergence + a_jump +
                                        a_pressure + a_dirichlet)        
        self.matrix_M = dolfin.assemble(dot(u, v)*dx)
        self.matrix_A = dolfin.assemble(a_diffusion)
        
        # These are not pre-assembled
        self.form_a_convection = a_convection
        self.form_L = L
    
    def assemble_lhs(self):
        data = self.simulation.data
        
        # Pre-assembled matrices
        Q = self.matrix_Q
        M = self.matrix_M
        A = self.matrix_A
        
        # Assemble convection matrix
        C = dolfin.assemble(self.form_a_convection)
        
        c1 = data['time_coeffs_py'][0]
        dt = 0.1
        Gr = 1.0 #float(data['rho_star'])
        Gm = 0.01 # float(data['nu_star'])/float(data['rho_star'])
        
        return Gr*c1/dt*M + Gm*A + Q + C
    
    def assemble_rhs(self):
        return dolfin.assemble(self.form_L)

# Default version
CoupledEquations = CoupledEquationsSlow#Fast

