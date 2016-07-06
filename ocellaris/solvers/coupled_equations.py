# encoding: utf8
import dolfin
from dolfin import dx, div, grad, dot, jump, avg, ds, dS, Constant
from . import UPWIND
from ..solver_parts import define_penalty


class CoupledEquations(object):
    def __init__(self, simulation, timestepping_method, flux_type, use_stress_divergence_form,
                 use_grad_p_form, use_grad_q_form, use_lagrange_multiplicator, 
                 pressure_continuity_factor, velocity_continuity_factor_D12,
                 include_hydrostatic_pressure, incompressibility_flux_type):
        """
        This class assembles the coupled Navier-Stokes equations, both CG and DG
        
        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
        self.timestepping_method = timestepping_method
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.use_grad_q_form = use_grad_q_form
        self.flux_type = flux_type
        self.use_lagrange_multiplicator = use_lagrange_multiplicator
        self.pressure_continuity_factor =  pressure_continuity_factor
        self.velocity_continuity_factor_D12 = velocity_continuity_factor_D12
        self.include_hydrostatic_pressure = include_hydrostatic_pressure
        self.incompressibility_flux_type = incompressibility_flux_type

        assert self.incompressibility_flux_type in ('central', 'upwind')
        
        # Discontinuous or continuous elements
        Vu_family = simulation.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_coupled_equation()
        
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
        self.simulation.log.info('DG SIP penalty:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
        if self.velocity_continuity_factor_D12 is not None:
            D12 = Constant([self.velocity_continuity_factor_D12]*self.simulation.ndim)
        else:
            D12 = Constant([0, 0])
        
        if self.pressure_continuity_factor != 0:
            h = self.simulation.data['h']
            h = Constant(1.0)
            D11 = avg(h/nu)*Constant(self.pressure_continuity_factor)
        else:
            D11 = None
        
        return Constant(penalty_dS), Constant(penalty_ds), D11, D12
    
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mpm = sim.multi_phase_model
        mesh = sim.data['mesh']
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
        
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties at t^{n}, t^{n-1} and t^{n-2}
        rho = mpm.get_density(0)
        rho_p = mpm.get_density(-1)
        rho_pp = mpm.get_density(-2)
        nu = mpm.get_laminar_kinematic_viscosity(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)
        rho_star = (rho + rho_p)/2
        
        # Hydrostatic pressure correction
        if self.include_hydrostatic_pressure:
            p += sim.data['p_hydrostatic']
        
        # Start building the coupled equations
        eq = 0
        
        # ALE mesh velocities
        u_mesh = dolfin.Constant([0]*sim.ndim)
        if sim.mesh_morpher.active:
            u_mesh = sim.data['u_mesh']
            
            # Modification of the convective velocity
            #u_conv -= u_mesh
            eq -= dot(div(rho*dolfin.outer(u, u_mesh)), v)*dx
            
            # Divergence of u should balance expansion/contraction of the cell K
            # ∇⋅u = -∂x/∂t       (See below for definition of the ∇⋅u term)
            cvol_new = dolfin.CellVolume(mesh)
            cvol_old = sim.data['cvolp']  
            eq += (cvol_new - cvol_old)/dt*q*dx
        
        if self.vel_is_discontinuous:
            penalty_dS, penalty_ds, D11, D12 = self.calculate_penalties(nu)
            
            # Upwind and downwind velocities
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        if self.use_lagrange_multiplicator:
            lm_trial = uc[ndim+1]
            lm_test = vc[ndim+1]
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
                eq += (rho_star*c1*u[d] + rho_p*c2*up + rho_pp*c3*upp)/dt*v[d]*dx
                
                # Convection
                # ∇⋅(ρ u ⊗ u_conv)
                eq += div(rho*u[d]*u_conv)*v[d]*dx
                
                if sim.mesh_morpher.active:
                    ud = up
                    um = -u_mesh
                    eq += div(rho*ud*um)*v[d]*dx
                
                # Diffusion
                # -∇⋅μ(∇u)
                eq += mu*dot(grad(u[d]), grad(v[d]))*dx
                
                # -∇⋅μ(∇u)^T
                if self.use_stress_divergence_form:
                    eq += mu*dot(u.dx(d), grad(v[d]))*dx
                
                # Pressure
                # ∇p
                eq -= v[d].dx(d)*p*dx
                
                # Body force (gravity)
                # ρ g
                eq -= rho*g[d]*v[d]*dx
                
                # Other sources
                for f in sim.data['momentum_sources']:
                    eq -= f[d]*v[d]*dx
                
                # Neumann boundary conditions
                neumann_bcs_pressure = sim.data['neumann_bcs'].get('p', [])
                for nbc in neumann_bcs_pressure:
                    eq += p*v[d]*n[d]*nbc.ds()
                
            else:
                # Weak form of the Navier-Stokes eq. with discontinuous elements
                assert self.flux_type == UPWIND
                
                # Divergence free criterion
                # ∇⋅u = 0
                if self.incompressibility_flux_type == 'central':
                    u_hat_p = avg(u[d])
                elif self.incompressibility_flux_type == 'upwind':
                    assert self.use_grad_q_form, 'Upwind only implemented for grad_q_form'
                    switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), 1.0, 0.0)
                    u_hat_p = switch*u[d]('+') + (1 - switch)*u[d]('-')
                
                if self.use_grad_q_form:
                    eq -= u[d]*q.dx(d)*dx
                    eq += (u_hat_p + D12[d]*jump(u, n))*jump(q)*n[d]('+')*dS
                else:
                    eq += q*u[d].dx(d)*dx
                    eq -= (avg(q) - dot(D12, jump(q, n)))*jump(u[d])*n[d]('+')*dS
                
                # Time derivative
                # ∂(ρu)/∂t
                eq += (rho_star*c1*u[d] + rho_p*c2*up + rho_pp*c3*upp)/dt*v[d]*dx
                
                # Convection:
                # -w⋅∇(ρu)
                flux_nU = rho*u[d]*w_nU
                flux = jump(flux_nU)
                eq -= rho*u[d]*div(v[d]*u_conv)*dx
                eq += flux*jump(v[d])*dS
                
                # ALE terms
                if sim.mesh_morpher.active:
                    ud = u[d]
                    um = -u_mesh
                    u_mesh_nU = (dot(um, n) + abs(dot(um, n)))/2.0
                    flux_mesh_nU = rho*ud*u_mesh_nU
                    flux_mesh = jump(flux_mesh_nU)
                    eq -= rho*ud*div(v[d]*um)*dx
                    eq += flux_mesh*jump(v[d])*dS
                
                # Diffusion:
                # -∇⋅∇u
                eq += mu*dot(grad(u[d]), grad(v[d]))*dx
                
                # Symmetric Interior Penalty method for -∇⋅μ∇u
                eq -= avg(mu)*dot(n('+'), avg(grad(u[d])))*jump(v[d])*dS
                eq -= avg(mu)*dot(n('+'), avg(grad(v[d])))*jump(u[d])*dS
                
                # Symmetric Interior Penalty coercivity term
                eq += penalty_dS*jump(u[d])*jump(v[d])*dS
                
                # -∇⋅μ(∇u)^T
                if self.use_stress_divergence_form:
                    eq += mu*dot(u.dx(d), grad(v[d]))*dx
                    eq -= avg(mu)*dot(n('+'), avg(u.dx(d)))*jump(v[d])*dS
                    eq -= avg(mu)*dot(n('+'), avg(v.dx(d)))*jump(u[d])*dS
                
                # Pressure
                # ∇p
                if self.use_grad_p_form:
                    eq += v[d]*p.dx(d)*dx
                    eq -= (avg(v[d]) + D12[d]*jump(v, n))*jump(p)*n[d]('+')*dS
                else:
                    eq -= p*v[d].dx(d)*dx
                    eq += (avg(p) - dot(D12, jump(p, n)))*jump(v[d])*n[d]('+')*dS
                
                # Pressure continuity stabilization. Needed for equal order discretization
                if D11 is not None:
                    eq += D11*dot(jump(p, n), jump(q, n))*dS
                
                # Body force (gravity)
                # ρ g
                eq -= rho*g[d]*v[d]*dx
                
                # Other sources
                for f in sim.data['momentum_sources']:
                    eq -= f[d]*v[d]*dx
                
                # Dirichlet boundary
                dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
                for dbc in dirichlet_bcs:
                    u_bc = dbc.func()
                    
                    # Divergence free criterion
                    if self.use_grad_q_form:
                        eq += q*u_bc*n[d]*dbc.ds()
                    else:
                        eq -= q*u[d]*n[d]*dbc.ds()
                        eq += q*u_bc*n[d]*dbc.ds()
                    
                    # Convection
                    eq += rho*u[d]*w_nU*v[d]*dbc.ds()
                    eq += rho*u_bc*w_nD*v[d]*dbc.ds()
                    
                    # SIPG for -∇⋅μ∇u
                    eq -= mu*dot(n, grad(u[d]))*v[d]*dbc.ds()
                    eq -= mu*dot(n, grad(v[d]))*u[d]*dbc.ds()
                    eq += mu*dot(n, grad(v[d]))*u_bc*dbc.ds()
                    
                    # Weak Dirichlet
                    eq += penalty_ds*(u[d] - u_bc)*v[d]*dbc.ds()
                    
                    # Pressure
                    if not self.use_grad_p_form:
                        eq += p*v[d]*n[d]*dbc.ds()
                
                # Neumann boundary
                neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
                for nbc in neumann_bcs:
                    # Divergence free criterion
                    if self.use_grad_q_form:
                        eq += q*u[d]*n[d]*nbc.ds()
                    else:
                        eq -= q*u[d]*n[d]*nbc.ds()
                    
                    # Convection
                    eq += rho*u[d]*w_nU*v[d]*nbc.ds()
                    
                    # Diffusion
                    eq -= mu*nbc.func()*v[d]*nbc.ds()
                    
                    # Pressure
                    if not self.use_grad_p_form:
                        eq += p*v[d]*n[d]*nbc.ds()
                
        a, L = dolfin.system(eq)
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


class CoupledEquationsScaledPressure(CoupledEquations):
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        Here we use the form with nu and 1/rho*grad(p)
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
        
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties
        rho = sim.data['rho_star']
        nu = sim.data['nu_star']
        phi = p/rho
        
        assert not self.use_stress_divergence_form
          
        if self.vel_is_discontinuous:
            penalty_dS, penalty_ds, D11, D12 = self.calculate_penalties(nu)
            
            # Upwind and downwind velocitues
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        if self.use_lagrange_multiplicator:
            lm_trial = uc[ndim+1]
            lm_test = vc[ndim+1]
            eq = (p*lm_test + q*lm_trial)*dx
        else:
            eq = 0
        
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
                eq += (c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
                
                # Convection
                # ∇⋅(u ⊗ w)
                eq += div(u[d]*u_conv)*v[d]*dx
                
                # Diffusion
                # -∇⋅ν∇u
                eq += nu*dot(grad(u[d]), grad(v[d]))*dx
                
                # Pressure
                # ∇ϕ
                eq -= phi*v[d].dx(d)*dx
                
                # Body force (gravity)
                # g
                eq -= g[d]*v[d]*dx
                
                # Neumann boundary conditions
                neumann_bcs_pressure = sim.data['neumann_bcs'].get('p', [])
                for nbc in neumann_bcs_pressure:
                    eq += phi*v[d]*n[d]*nbc.ds()
                
            else:
                # Weak form of the Navier-Stokes eq. with discontinuous elements
                assert self.flux_type == UPWIND
                
                # Divergence free criterion
                # ∇⋅u = 0
                if self.use_grad_p_form:
                    eq -= u[d]*q.dx(d)*dx
                    eq += (avg(u[d]) + D12[d]*jump(u, n))*jump(q)*n[d]('+')*dS
                else:
                    eq += q*u[d].dx(d)*dx
                    eq -= (avg(q) - dot(D12, jump(q, n)))*jump(u[d])*n[d]('+')*dS
                
                # Time derivative
                # ∂u/∂t
                eq += (c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
                
                # Convection:
                # -w⋅∇u    
                flux_nU = u[d]*w_nU
                flux = jump(flux_nU)
                eq -= u[d]*div(v[d]*u_conv)*dx
                eq += flux*jump(v[d])*dS
                eq += flux_nU*v[d]*ds
                
                # Diffusion:
                # -∇⋅ν∇u
                eq += nu*dot(grad(u[d]), grad(v[d]))*dx
                
                # Symmetric Interior Penalty method for -∇⋅μ∇u
                eq -= avg(nu)*dot(n('+'), avg(grad(u[d])))*jump(v[d])*dS
                eq -= avg(nu)*dot(n('+'), avg(grad(v[d])))*jump(u[d])*dS
                
                # Symmetric Interior Penalty coercivity term
                eq += penalty_dS*jump(u[d])*jump(v[d])*dS
                
                # Pressure
                # ∇ϕ
                if self.use_grad_p_form:
                    eq += v[d]*phi.dx(d)*dx
                    eq -= (avg(v[d]) + D12[d]*jump(v, n))*jump(phi)*n[d]('+')*dS
                else:
                    eq -= phi*v[d].dx(d)*dx
                    eq += (avg(phi) - dot(D12, jump(p, n)))*jump(v[d])*n[d]('+')*dS
                
                # Pressure continuity stabilization. Needed for equal order discretization    
                if D11 is not None:
                    eq += D11*dot(jump(p, n), jump(q, n))*dS
                
                # Body force (gravity)
                # ρ g
                eq -= g[d]*v[d]*dx
                
                # Dirichlet boundary
                dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
                for dbc in dirichlet_bcs:
                    u_bc = dbc.func()
                    
                    # Divergence free criterion
                    if self.use_grad_p_form:
                        eq += q*u_bc*n[d]*dbc.ds()
                    else:
                        eq -= q*u[d]*n[d]*dbc.ds()
                        eq += q*u_bc*n[d]*dbc.ds()
                    
                    # Convection
                    eq += w_nD*dbc.func()*v[d]*dbc.ds()
                    
                    # SIPG for -∇⋅μ∇u
                    eq -= nu*dot(n, grad(u[d]))*v[d]*dbc.ds()
                    eq -= nu*dot(n, grad(v[d]))*u[d]*dbc.ds()
                    eq += nu*dot(n, grad(v[d]))*u_bc*dbc.ds()
                    
                    # Weak Dirichlet
                    eq += penalty_ds*(u[d] - u_bc)*v[d]*dbc.ds()
                    
                    # Pressure
                    eq += phi*v[d]*n[d]*dbc.ds()
        
        a, L = dolfin.system(eq)
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


class CoupledEquationsPreassembled(CoupledEquations):
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        
        This implementation tries to avoid re-assembly of matrices as much as possible
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
            penalty_dS, penalty_ds = self.calculate_penalties()
            
            # Upwind and downwind velocitues
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        a_lm = (p*lm_test + q*lm_trial)*dx

        # Forms to be assembled for each dimension
        a_divergence = 0
        a_convection = 0
        a_diffusion = 0
        a_jump = 0
        a_pressure = 0
        a_dirichlet = 0
        L = 0
        
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


EQUATION_SUBTYPES = {
    'Conservative': CoupledEquations,
    'ScaledPressure': CoupledEquationsScaledPressure,
    #'ConservativePreassembled': CoupledEquationsPreassembled,
}