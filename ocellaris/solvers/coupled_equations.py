# encoding: utf8
from __future__ import division
import dolfin
from dolfin import dx, div, grad, dot, jump, avg, dS, Constant
from . import UPWIND
from ..solver_parts import define_penalty


class CoupledEquations(object):
    def __init__(self, simulation, flux_type, use_stress_divergence_form,
                 use_grad_p_form, use_grad_q_form, use_lagrange_multiplicator, 
                 pressure_continuity_factor, velocity_continuity_factor_D12,
                 include_hydrostatic_pressure, incompressibility_flux_type):
        """
        This class assembles the coupled Navier-Stokes equations, both CG and DG
        
        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
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
        self.simulation.log.info('    DG SIP penalty:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
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
        
        # Fluid properties
        rho = mpm.get_density(0)
        nu = mpm.get_laminar_kinematic_viscosity(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)
        
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
                eq += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
                
                # Convection
                # ∇⋅(ρ u ⊗ u_conv)
                eq += rho*dot(u_conv, grad(u[d]))*v[d]*dx
                
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
                    
                # Outlet boundary
                for obc in sim.data['outlet_bcs']:     
                    # Diffusion
                    mu_dudn = p*n[d]
                    eq -= mu_dudn*v[d]*obc.ds()
                       
                    # Pressure
                    p_ = mu*dot(dot(grad(u), n), n)
                    eq += p_*v[d]*n[d]*obc.ds()
            
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
                eq += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
                
                # Convection:
                # -w⋅∇(ρu)
                flux_nU = u[d]*w_nU
                flux = jump(flux_nU)
                eq -= u[d]*dot(grad(rho*v[d]), u_conv)*dx
                eq += flux*jump(rho*v[d])*dS
                
                # Stabilizing term when w is not divergence free
                eq += 1/2*div(u_conv)*u[d]*v[d]*dx
                
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
                
                # Outlet boundary
                for obc in sim.data['outlet_bcs']:
                    # Divergence free criterion
                    if self.use_grad_q_form:
                        eq += q*u[d]*n[d]*obc.ds()
                    else:
                        eq -= q*u[d]*n[d]*obc.ds()
                    
                    # Convection
                    eq += rho*u[d]*w_nU*v[d]*obc.ds()
                    
                    # Diffusion
                    mu_dudn = p*n[d]
                    eq -= mu_dudn*v[d]*obc.ds()
                    
                    # Pressure
                    if not self.use_grad_p_form:
                        p_ = mu*dot(dot(grad(u), n), n)
                        eq += p_*n[d]*v[d]*obc.ds()
                
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


class CoupledEquationsPoissonPressure(object):
    def __init__(self, simulation, flux_type, use_stress_divergence_form,
                 use_grad_p_form, use_grad_q_form, use_lagrange_multiplicator, 
                 pressure_continuity_factor, velocity_continuity_factor_D12,
                 include_hydrostatic_pressure, incompressibility_flux_type):
        """
        This class assembles the coupled Navier-Stokes equations with a Poisson equation
        for the pressure
        
            ρ∂u/∂t + ρw⋅∇u = ∇⋅μ∇u -∇p + ρg
                          -∇⋅∇p = ρ∇w:∇u^T 
        
        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
        self.include_hydrostatic_pressure = include_hydrostatic_pressure
        self.use_grad_p_form = use_grad_p_form
        
        assert not use_stress_divergence_form
        assert not use_lagrange_multiplicator
        assert not simulation.mesh_morpher.active
        assert incompressibility_flux_type == 'central'
        assert simulation.data['Vu'].ufl_element().family() == 'Discontinuous Lagrange'
        
        # Create UFL forms
        self.define_coupled_equation()
    
    def calculate_penalties(self):
        """
        Calculate SIPG penalty
        """
        mpm = self.simulation.multi_phase_model
        mesh = self.simulation.data['mesh']
        
        mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
        Pu = self.simulation.data['Vu'].ufl_element().degree()
        Pp = self.simulation.data['Vp'].ufl_element().degree()
        penalty_u_dS = define_penalty(mesh, Pu, mu_min, mu_max, boost_factor=3, exponent=1.0)
        penalty_p_dS = define_penalty(mesh, Pp, 1.0, 1.0, boost_factor=3, exponent=1.0)
        penalty_u_ds = penalty_u_dS*2
        penalty_p_ds = penalty_p_dS*2
        self.simulation.log.info('    DG SIP penalty u:  dS %.1f  ds %.1f' % (penalty_u_dS, penalty_u_ds))
        self.simulation.log.info('    DG SIP penalty p:  dS %.1f  ds %.1f' % (penalty_p_dS, penalty_p_ds))
        
        return Constant(penalty_u_dS), Constant(penalty_u_ds), Constant(penalty_p_dS), Constant(penalty_p_ds)
    
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
        
        # Fluid properties
        rho = mpm.get_density(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)
        
        # Hydrostatic pressure correction
        if self.include_hydrostatic_pressure:
            p += sim.data['p_hydrostatic']
        
        # Start building the coupled equations
        eq = 0
        
        penalty_u_dS, penalty_u_ds, penalty_p_dS, _penalty_p_ds = self.calculate_penalties()
        
        # Upwind and downwind velocities
        w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
        w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
    
        # Momentum equations
        for d in range(sim.ndim):
            up = sim.data['up%d' % d]
            upp = sim.data['upp%d' % d]
            
            # Time derivative
            # ρ∂u/∂t
            eq += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
            
            # Convection:
            # ρw⋅∇u
            flux_nU = u[d]*w_nU
            flux = jump(flux_nU)
            eq -= u[d]*dot(grad(rho*v[d]), u_conv)*dx
            eq += flux*jump(rho*v[d])*dS
            
            # Stabilizing term when w is not divergence free
            eq += 1/2*div(u_conv)*u[d]*v[d]*dx
            
            # Diffusion:
            # -∇⋅μ∇u
            eq += mu*dot(grad(u[d]), grad(v[d]))*dx
            
            # Symmetric Interior Penalty method for -∇⋅μ∇u
            eq -= avg(mu)*dot(n('+'), avg(grad(u[d])))*jump(v[d])*dS
            eq -= avg(mu)*dot(n('+'), avg(grad(v[d])))*jump(u[d])*dS
            
            # Symmetric Interior Penalty coercivity term
            eq += penalty_u_dS*jump(u[d])*jump(v[d])*dS
            
            # Pressure
            # ∇p
            if self.use_grad_p_form:
                eq += v[d]*p.dx(d)*dx
                eq -= avg(v[d])*jump(p)*n[d]('+')*dS
            else:
                eq -= p*v[d].dx(d)*dx
                eq += avg(p)*jump(v[d])*n[d]('+')*dS
            
            # Body force (gravity)
            # ρ g
            eq -= rho*g[d]*v[d]*dx
            
            # Other sources
            for f in sim.data['momentum_sources']:
                eq -= f[d]*v[d]*dx
                
            # Pressure Poisson equation
            # -∇⋅∇p = ρ∇w:∇u^T
            eq += dot(grad(p), grad(q))*dx

            # Symmetric Interior Penalty method for -∇⋅∇p
            eq -= dot(n('+'), avg(grad(p)))*jump(q)*dS
            eq -= dot(n('+'), avg(grad(q)))*jump(p)*dS

            # Symmetric Interior Penalty coercivity term
            eq += penalty_p_dS*jump(p)*jump(q)*dS
            
            # -ρ∇w:∇u^T
            eq -= dot(avg(u_conv.dx(d)), n('+'))*avg(u[d])*jump(rho*q)*dS
            eq += dot(u_conv.dx(d), grad(rho*q))*u[d]*dx
                        
            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                
                # Convection
                eq += rho*u[d]*w_nU*v[d]*dbc.ds()
                eq += rho*u_bc*w_nD*v[d]*dbc.ds()
                
                # SIPG for -∇⋅μ∇u
                eq -= mu*dot(n, grad(u[d]))*v[d]*dbc.ds()
                eq -= mu*dot(n, grad(v[d]))*u[d]*dbc.ds()
                eq += mu*dot(n, grad(v[d]))*u_bc*dbc.ds()
                
                # Weak Dirichlet
                eq += penalty_u_ds*(u[d] - u_bc)*v[d]*dbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    eq += p*v[d]*n[d]*dbc.ds()
                    
                # Pressure Poisson equation
                #dpdn = dot(n, grad(p))
                if d == 0:
                    dpdn = dot(n, mu*div(grad(u))
                                  + rho*g
                                  - rho*dot(grad(u), u_conv)
                                  - rho*(c1*u + c2*sim.data['up'] + c3*sim.data['upp'])/dt) 
                    eq -= dpdn*q*dbc.ds()
                eq -= rho*dot(u_conv.dx(d), n)*u_bc*q*dbc.ds()
            
            # Neumann boundary
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
            for nbc in neumann_bcs:
                assert False
                # Convection
                eq += rho*u[d]*w_nU*v[d]*nbc.ds()
                
                # Diffusion
                eq -= mu*nbc.func()*v[d]*nbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    eq += p*v[d]*n[d]*nbc.ds()
                    
                # Pressure Poisson equation
                eq -= dot(n, grad(p))*q*dbc.ds()
                eq -= dot(u_conv.dx(d), n)*u[d]*q*dbc.ds()
            
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


EQUATION_SUBTYPES = {
    'Default': CoupledEquations,
    'PoissonP': CoupledEquationsPoissonPressure,
}