# encoding: utf8
from __future__ import division
import dolfin
from dolfin import dot, div, grad, avg, jump, dx, dS, Constant
from . import BDF, UPWIND
from ..solver_parts import define_penalty


class BaseEquation(object):
    # Will be shadowed by object properties after first assemble
    tensor_lhs = None
    tensor_rhs = None
    
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


class MomentumPredictionEquation(BaseEquation):
    def __init__(self, simulation, component, timestepping_method, flux_type,
                 use_stress_divergence_form, use_grad_p_form, include_hydrostatic_pressure):
        """
        This class assembles the momentum equation for one velocity component, both CG and DG 
        """
        self.simulation = simulation
        self.component = component
        self.timestepping_method = timestepping_method
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.flux_type = flux_type
        self.include_hydrostatic_pressure = include_hydrostatic_pressure
        
        assert self.timestepping_method == BDF
        
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
        self.simulation.log.info('    DG SIP penalty viscosity:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
        D12 = Constant([0, 0])
        
        return Constant(penalty_dS), Constant(penalty_ds), D12
    
    def define_momentum_equation(self):
        """
        Setup the momentum equation for one velocity component
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mpm = sim.multi_phase_model
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
        
        # ALE mesh velocities
        if sim.mesh_morpher.active:
            u_mesh = sim.data['u_mesh']
            u_conv -= u_mesh
        
        if self.include_hydrostatic_pressure:
            p += sim.data['p_hydrostatic']
        
        # Fluid properties
        rho = sim.data['rho']
        nu = mpm.get_laminar_kinematic_viscosity(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)
        
        # Include (∇u)^T term?
        assert not self.use_stress_divergence_form
        
        # Values at previous time steps
        up = sim.data['up%d' % self.component]
        upp = sim.data['upp%d' % self.component]
            
        if not self.vel_is_discontinuous:
            # Weak form of the Navier-Stokes eq. with continuous elements
            
            # Time derivative
            # ∂u/∂t
            a = rho*c1*u/dt*v*dx
            L = -(rho*c2*up + rho*c3*upp)/dt*v*dx
            
            # Convection
            # ρ∇⋅(u ⊗ u_conv)
            a += rho*div(u*u_conv)*v*dx
            
            # Diffusion
            # -∇⋅μ∇u
            a += mu*dot(grad(u), grad(v))*dx
            
            # Pressure
            # ∇p
            if self.use_grad_p_form:
                L -= p.dx(self.component)*v*dx
            else:
                L += v.dx(self.component)*p*dx
            
            # Body force (gravity)
            # ρ g
            L += rho*g[self.component]*v*dx
            
            # Neumann boundary
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % self.component, [])
            for nbc in neumann_bcs:
                L += mu*nbc.func()*v*nbc.ds()
                
                if not self.use_grad_p_form:
                    L -= p*v*ni*nbc.ds()
        
        else:
            # Weak form of the Navier-Stokes eq. with discontinuous elements
            assert self.flux_type == UPWIND
            
            # Upwind and downwind velocities
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
            
            # Penalties
            penalty_dS, penalty_ds, D12 = self.calculate_penalties(nu)
            
            # Time derivative
            # ρ ∂(u)/∂t
            a = rho*c1*u/dt*v*dx
            L = -rho*(c2*up + c3*upp)/dt*v*dx
            
            # Convection:
            # w⋅∇(ρu)    
            flux_nU = rho*u*w_nU
            flux = jump(flux_nU)
            a -= rho*u*div(v*u_conv)*dx
            a += flux*jump(v)*dS
            
            # Skew-symetric stabilizing term when w is not divergence free
            a += 1/2*div(u_conv)*u*v*dx
            
            # Diffusion:
            # -∇⋅∇u
            a += mu*dot(grad(u), grad(v))*dx
            
            # Symmetric Interior Penalty method for -∇⋅μ∇u
            a -= avg(mu)*dot(n('+'), avg(grad(u)))*jump(v)*dS
            a -= avg(mu)*dot(n('+'), avg(grad(v)))*jump(u)*dS
            
            # Symmetric Interior Penalty coercivity term
            a += penalty_dS*jump(u)*jump(v)*dS
            
            # Pressure
            # ∇p
            if self.use_grad_p_form:
                L -= v*p.dx(self.component)*dx
                L += (avg(v) + dot(D12, jump(v, n)))*jump(p)*ni('+')*dS
            else:
                L += p*v.dx(self.component)*dx
                L -= (avg(p) - dot(D12, jump(p, n)))*jump(v)*ni('+')*dS
            
            # Body force (gravity)
            # ρ g
            L += rho*g[self.component]*v*dx
            
            # Other sources
            for f in sim.data['momentum_sources']:
                L+= f[self.component]*v*dx
            
            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % self.component, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                
                # Convection
                a += rho*u*w_nU*v*dbc.ds()
                L -= rho*u_bc*w_nD*v*dbc.ds()
                
                # SIPG for -∇⋅μ∇u
                a -= mu*dot(n, grad(u))*v*dbc.ds()
                a -= mu*dot(n, grad(v))*u*dbc.ds()
                L -= mu*dot(n, grad(v))*u_bc*dbc.ds()
                
                # Weak Dirichlet
                a += penalty_ds*u*v*dbc.ds()
                L += penalty_ds*u_bc*v*dbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    L -= p*v*ni*dbc.ds()
            
            # Neumann boundary
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % self.component, [])
            for nbc in neumann_bcs:
                # Convection
                a += rho*u*w_nU*v*nbc.ds()
                
                # Diffusion
                L += mu*nbc.func()*v*nbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    L -= p*v*ni*nbc.ds()
            
            # Outlet boundary
            for obc in sim.data['outlet_bcs']:
                # Convection
                a += rho*u*w_nU*v*obc.ds()
                
                # Diffusion
                mu_dudn = p*ni
                L += mu_dudn*v*obc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    L -= p*v*ni*obc.ds()
        
        self.form_lhs = a
        self.form_rhs = L


class PressureCorrectionEquation(BaseEquation):
    def __init__(self, simulation, use_lagrange_multiplicator, incompressibility_flux_type):
        """
        This class assembles the pressure Poisson equation, both CG and DG 
        """
        self.simulation = simulation
        self.use_lagrange_multiplicator = use_lagrange_multiplicator
        self.incompressibility_flux_type = incompressibility_flux_type
        
        # Discontinuous or continuous elements
        Vp_family = simulation.data['Vp'].ufl_element().family()
        self.pressure_is_discontinuous = (Vp_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_pressure_equation()
        
        self.matrix = None
        self.pressure_null_space = None
        
    def solve(self, solver, assemble_A, remove_null_space):
        """
        Solve the pressure equation
        """
        p = self.simulation.data['p']
        if assemble_A:
            self.matrix = self.assemble_lhs()
        
        # The equation system to solve
        A = self.matrix
        b = self.assemble_rhs()
        
        # Apply strong boundary conditions
        if not self.pressure_is_discontinuous:
            dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
            for dbc in dirichlet_bcs:
                dbc.apply(A, b)
        
        # Inform PETSc about the null space
        if remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p.vector())
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                
                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])
            
            # Make sure the null space is set on the matrix
            dolfin.as_backend_type(A).set_nullspace(self.pressure_null_space)
            
            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(b)
        
        # Solve for the new pressure
        niters_p = solver.solve(A, p.vector(), b)
        
        return niters_p
    
    def calculate_penalties(self):
        """
        Calculate SIPG penalty
        """
        mesh = self.simulation.data['mesh']
        P = self.simulation.data['Vp'].ufl_element().degree()
        rho_min, rho_max = self.simulation.multi_phase_model.get_density_range()
        k_min = 1.0/rho_max
        k_max = 1.0/rho_min
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
        u_star = sim.data['u_star']
        
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
        rho = sim.data['rho']
        
        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        assert not self.use_lagrange_multiplicator, 'NOT IMPLEMENTED YET'
        
        if not self.pressure_is_discontinuous:
            # Weak form of the Poisson eq. with continuous elements
            # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^* 
            a = dot(grad(p), grad(q))*dx
            L = dot(grad(p_star), grad(q))*dx
            L -= c1/dt*rho*div(u_star)*q*dx
            
            # Neumann boundary conditions on p and p_star cancel
        
        else:
            # Weak form of the Poisson eq. with discontinuous elements
            # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^*
            K = 1.0/rho
            a = K*dot(grad(p), grad(q))*dx
            L = K*dot(grad(p_star), grad(q))*dx

            # RHS, ∇⋅u^*
            if self.incompressibility_flux_type == 'central':
                u_flux = avg(u_star)
            elif self.incompressibility_flux_type == 'upwind':
                switch = dolfin.conditional(dolfin.gt(abs(dot(u_star, n))('+'), 0.0), 1.0, 0.0)
                u_flux = switch*u_star('+') + (1 - switch)*u_star('-')
            L += c1/dt*dot(u_star, grad(q))*dx
            L -= c1/dt*dot(u_flux, n('+'))*jump(q)*dS
            
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
            
            # Neumann boundary conditions
            neumann_bcs = sim.data['neumann_bcs'].get('p', [])
            for nbc in neumann_bcs:
                # Neumann boundary conditions on p and p_star cancel
                #L += (nbc.func() - dot(n, grad(p_star)))*q*nbc.ds()
                pass
            
            # Use boundary conditions for the velocity for the
            # term from integration by parts of div(u_star)
            for d in range(sim.ndim):
                dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
                neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
                for dbc in dirichlet_bcs:
                    u_bc = dbc.func()
                    L -= c1/dt*u_bc*n[d]*q*dbc.ds()
                for nbc in neumann_bcs:
                    L -= c1/dt*u_star[d]*n[d]*q*nbc.ds()
        
        # ALE mesh velocities
        if sim.mesh_morpher.active:
            cvol_new = dolfin.CellVolume(mesh)
            cvol_old = sim.data['cvolp']
            
            # Divergence of u should balance expansion/contraction of the cell K
            # ∇⋅u = -∂K/∂t       (See below for definition of the ∇⋅u term)  
            L -= (cvol_new - cvol_old)/dt*q*dx
        
        self.form_lhs = a
        self.form_rhs = L


class LDGPressureCorrectionEquation(BaseEquation):
    def __init__(self, simulation, use_lagrange_multiplicator, incompressibility_flux_type):
        """
        This class assembles the pressure Poisson equation, DG LDG method 
        """
        self.simulation = simulation
        assert not use_lagrange_multiplicator
        assert incompressibility_flux_type == 'central'
        assert simulation.data['Vp'].ufl_element().family() == 'Discontinuous Lagrange'
        assert not simulation.mesh_morpher.active
        
        mesh = simulation.data['mesh']
        Vp = simulation.data['Vp']
        P = Vp.ufl_element().degree()
        ep = Vp.ufl_element()
        es = dolfin.VectorElement('DG', mesh.ufl_cell(), P)
        em = dolfin.MixedElement([ep, es])
        self.Vcoupled = dolfin.FunctionSpace(mesh, em)
        
        # Create UFL forms
        self.define_pressure_equation()
        
        self.assigner = dolfin.FunctionAssigner(simulation.data['Vp'], self.Vcoupled.sub(0))
        self.coupled_func = dolfin.Function(self.Vcoupled)
        self.matrix = None
        self.pressure_null_space = None
    
    def solve(self, solver, assemble_A, remove_null_space):
        """
        Solve the pressure equation
        """
        p = self.simulation.data['p']
        if assemble_A:
            self.matrix = self.assemble_lhs()
        
        # The equation system to solve
        A = self.matrix
        b = self.assemble_rhs()
        
        # Inform PETSc about the null space
        if remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p.vector())
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                
                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])
                
                # Create null space vector in Vp Space
                null_func = dolfin.Function(self.simulation.data['Vp'])
                null_vec = null_func.vector()
                null_vec[:] = 1
                null_vec *= 1 / null_vec.norm("l2")
                
                # Convert null space vector to coupled space
                null_func2 = dolfin.Function(self.Vcoupled)
                fa = dolfin.FunctionAssigner(self.Vcoupled.sub(0), self.simulation.data['Vp'])
                fa.assign(null_func2.sub(0), null_func)
                
                # Create the null space basis
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_func2.vector()])
            
            # Make sure the null space is set on the matrix
            dolfin.as_backend_type(A).set_nullspace(self.pressure_null_space)
            
            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(b)
        
        # Solve for the new pressure
        niters_p = solver.solve(A, self.coupled_func.vector(), b)
        self.assigner.assign(p, self.coupled_func.sub(0))
        
        return niters_p
    
    def define_pressure_equation(self):
        """
        Setup the pressure Poisson equation
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        p_star = sim.data['p']
        u_star = sim.data['u_star']
        
        # Trial and test functions
        p, sigma = dolfin.TrialFunctions(self.Vcoupled)
        q, tau = dolfin.TestFunctions(self.Vcoupled)
        
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties
        mpm = sim.multi_phase_model
        mu = mpm.get_laminar_dynamic_viscosity(0)
        rho = sim.data['rho']
        
        # LDG interior facet fluxes
        C11 = Constant(1.0)
        C12 = 0.5*(0.5*n('+') + 0.5*n('-'))
        p_hat_dS = avg(p) - dot(C12, jump(p, n))
        sigma_hat_dS = avg(sigma) - C11*jump(p, n) + C12*jump(sigma, n)
        
        # Weak form of the Poisson eq. with discontinuous elements
        # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^*
        K = 1.0/rho
        
        # Equation 1
        eq = dot(sigma, tau)*dx
        eq += p*div(K*tau)*dx
        eq -= p_hat_dS*jump(K*tau, n)*dS
        
        # Equation 2
        eq += dot(sigma, grad(q))*dx
        eq -= dot(sigma_hat_dS, jump(q, n))*dS
        
        # P star
        eq -= K*dot(grad(p_star), grad(q))*dx
        eq += dot(n('+'), avg(K*grad(p_star)))*jump(q)*dS
        
        # RHS, ∇⋅u^*
        u_flux = avg(u_star)
        eq -= c1/dt*dot(u_star, grad(q))*dx
        eq += c1/dt*dot(u_flux, n('+'))*jump(q)*dS
        
        # Collect Dirichlet and outlet boundary values
        dirichlet_vals_and_ds = []
        for dbc in sim.data['dirichlet_bcs'].get('p', []):
            dirichlet_vals_and_ds.append((dbc.func(), dbc.ds()))
        for obc in sim.data['outlet_bcs']:
            p_ = mu*dot(dot(grad(u_star), n), n)
            dirichlet_vals_and_ds.append((p_, obc.ds()))
        
        # Dirichlet boundary
        for p_bc, dds in dirichlet_vals_and_ds:
            # LDG exterior facet fluxes
            p_hat_ds = p_bc
            sigma_hat_ds = sigma - C11*(p - p_bc)*n
            
            # Equation 1
            eq -= p_hat_ds*K*dot(tau, n)*dds
            
            # Equation 2
            eq -= dot(sigma_hat_ds, n)*q*dds
            
            # P star
            eq += dot(n, K*grad(p_star))*q*dds
            
        neumann_bcs = sim.data['neumann_bcs'].get('p', [])
        for nbc in neumann_bcs:
            # LDG exterior facet fluxes
            p_hat_ds = p
            #sigma_hat_ds = nbc.func()*n
            
            # Equation 1
            eq -= p_hat_ds*K*dot(tau, n)*nbc.ds()
            
            # Equation 2
            #eq -= dot(sigma_hat_ds, n)*q*nbc.ds()
            
            # P star
            #eq += dot(n, K*grad(p_star))*q*nbc.ds()
        
        # Use boundary conditions for the velocity for the
        # term from integration by parts of div(u_star)
        for d in range(sim.ndim):
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                eq += c1/dt*u_bc*n[d]*q*dbc.ds()
            for nbc in neumann_bcs:
                eq += c1/dt*u_star[d]*n[d]*q*nbc.ds()
        
        a, L = dolfin.system(eq)
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
        rho = sim.data['rho']
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        
        Vu = sim.data['Vu']
        us = sim.data['u_star%d' % self.component]
        p_hat = sim.data['p_hat']
        u = dolfin.TrialFunction(Vu)
        v = dolfin.TestFunction(Vu)
        
        self.form_lhs = u*v*dx
        self.form_rhs = us*v*dx - dt/(c1*rho)*p_hat.dx(self.component)*v*dx


EQUATION_SUBTYPES = {
    'Default': (MomentumPredictionEquation, PressureCorrectionEquation, VelocityUpdateEquation),
    'LDGforPressure': (MomentumPredictionEquation, LDGPressureCorrectionEquation, VelocityUpdateEquation),
}
