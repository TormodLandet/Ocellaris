# encoding: utf8
from __future__ import division
import dolfin
from dolfin import dx, div, grad, dot, jump, avg, dS
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error, timeit, linear_solver_from_input
from . import Solver, register_solver, BDF, BLENDED, UPWIND, LOCAL_LAX_FRIEDRICH
from .ipcs_equations import define_penalty


# Default values, can be changed in the input file
SOLVER = 'petsc'
LU_PARAMETERS = {}

# Implemented timestepping methods
TIMESTEPPING_METHODS = (BDF,)


@register_solver('Coupled')
class SolverCoupled(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the pressure-velocity splitting
        scheme IPCS (Incremental Pressure Correction Scheme)
        """
        self.simulation = sim = simulation
        self.create_functions()
        self.read_input()
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        self.is_first_timestep = True
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        self.is_single_phase = isinstance(sim.data['rho'], dolfin.Constant)
        
        # Get the BCs for the coupled function space
        self.dirichlet_bcs = self.coupled_boundary_conditions()
        
        # Create equation
        self.define_coupled_equation()
        
        # Store number of iterations
        self.niters = None
    
    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation
        
        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        Vl = dolfin.FunctionSpace(sim.data['mesh'], "R", 0)
        Vu_family = Vu.ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
            
        # Create coupled mixed function space and mixed function to hold results
        if sim.ndim == 2:
            self.Vcoupled = Vu * Vu * Vp * Vl
            self.subspaces = [self.Vcoupled.sub(0).sub(0).sub(0),
                              self.Vcoupled.sub(0).sub(0).sub(1),
                              self.Vcoupled.sub(0).sub(1),
                              self.Vcoupled.sub(1)]
            self.subspace_names = 'u0 u1 p l'.split()
        
        else:
            self.Vcoupled = Vu * Vu * Vu * Vp * Vl  
            self.subspaces = [self.Vcoupled.sub(0).sub(0).sub(0),
                              self.Vcoupled.sub(0).sub(0).sub(1),
                              self.Vcoupled.sub(0).sub(0).sub(2),
                              self.Vcoupled.sub(0).sub(0).sub(1),
                              self.Vcoupled.sub(1)]
            self.subspace_names = 'u0 u1 u2 p l'.split()
        
        # Create segregated velocity functions on component and vector form
        u_list, up_list, upp_list, u_conv = [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            u_list.append(u)
            up_list.append(up)
            upp_list.append(upp)
            u_conv.append(uc)
        sim.data['u'] = dolfin.as_vector(u_list)
        sim.data['up'] = dolfin.as_vector(up_list)
        sim.data['upp'] = dolfin.as_vector(upp_list)
        sim.data['u_conv'] = dolfin.as_vector(u_conv)
        sim.data['coupled'] = self.coupled_func = dolfin.Function(self.Vcoupled)
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Solver for the coupled system
        self.coupled_solver = linear_solver_from_input(sim, 'solver/coupled', 'lu', 
                                                       None, SOLVER, LU_PARAMETERS)
        
        # Give warning if using iterative solver
        if isinstance(self.coupled_solver, dolfin.KrylovSolver):
            sim.log.warning('WARNING: Using a Krylov solver for the coupled NS equations is not a good idea')
        else:
            self.coupled_solver.parameters['same_nonzero_pattern'] = True
        
        # Get flux type for the velocity
        self.flux_type = sim.input.get_value('convection/u/flux_type', BLENDED, 'string')
        if self.vel_is_discontinuous and self.flux_type == BLENDED:
            # Get convection blending scheme
            conv_schemes = []
            conv_scheme_name = sim.input.get_value('convection/u/convection_scheme', 'Upwind', 'string')
            for d in range(sim.ndim):
                conv_scheme = get_convection_scheme(conv_scheme_name)(sim, 'u%d' % d)
                conv_schemes.append(conv_scheme)
            self.convection_schemes = conv_schemes
        
        # Coefficients for u, up and upp
        self.timestepping_method = sim.input.get_value('solver/timestepping_method', BDF, 'string')
        if not self.timestepping_method in TIMESTEPPING_METHODS:
            available_methods = '\n'.join(' - %s' % m for m in TIMESTEPPING_METHODS)
            report_error('Unknown timestepping method', 
                         'Timestepping method %s not recognised, please use one of:\n%s' %
                         (self.timestepping_method, available_methods))
    
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        """
        uc = dolfin.TrialFunction(self.Vcoupled)
        vc = dolfin.TestFunction(self.Vcoupled)
        
        # Unpack the coupled trial and test functions
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
        
        sim = self.simulation
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        u_conv = sim.data['u_conv']
        g = sim.data['g']
        rho = sim.data['rho']
        nu = sim.data['nu']
        mu = rho*nu
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Include (∇u)^T term?
        sd = dolfin.Constant(1.0)
        
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
                eq += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
                
                # Convection
                # ∇⋅(ρ u ⊗ u_conv)
                eq += div(rho*u[d]*u_conv)*v[d]*dx
                
                # Diffusion
                # -∇⋅μ[(∇u) + (∇u)^T]
                eq += mu*dot(grad(u[d]), grad(v[d]))*dx
                eq += sd*mu*dot(u.dx(d), grad(v[d]))*dx
                
                # Pressure
                # ∇p
                eq -= v[d].dx(d)*p*dx
                
                # Body force (gravity)
                # ρ g
                eq -= rho*g[d]*v[d]*dx
            
            else:
                # Weak form of the Navier-Stokes eq. with discontinuous elements
                
                # Divergence free criterion
                # ∇⋅u = 0
                eq += u[d].dx(d)*q*dx
                eq += avg(q)*jump(u[d])*n[d]('+')*dS
                
                # Time derivative
                # ∂u/∂t
                eq += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
                
                # Define the convective flux
                # f = ρ u ⊗ u_conv
                def calc_flux(ui, uc):
                    if self.flux_type == BLENDED:
                        # Upstream and downstream normal velocities
                        flux_nU = rho*ui*(dot(uc, n) + abs(dot(uc, n)))/2
                        flux_nD = rho*ui*(dot(uc, n) - abs(dot(uc, n)))/2
                        
                        # Define the blended flux
                        # The blending factor beta is not DG, so beta('+') == beta('-')
                        b = self.convection_schemes[d].blending_function('+')
                        flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
                        flux = flux_nU('+') - flux_nU('-')
                    
                    elif self.flux_type == UPWIND:
                        # Pure upwind flux
                        flux_nU = rho*ui*(dot(uc, n) + abs(dot(uc, n)))/2
                        flux = flux_nU('+') - flux_nU('-')
                    
                    elif self.flux_type == LOCAL_LAX_FRIEDRICH:
                        # Local Lax-Friedrich flux
                        uflmax = lambda a, b: dolfin.conditional(dolfin.gt(a, b), a, b)
                        uflmag = lambda a: dolfin.sqrt(dolfin.inner(a, a))
                        C = 1.2*uflmax(uflmag(rho('+')*uc('+')), uflmag(rho('-')*uc('-')))
                        flux = avg(rho*ui*uc) + C/2*jump(ui)*n('+')
                        flux = dot(flux, n('+'))
                    return flux
                
                # Convection
                # ∇⋅(ρ u ⊗ u_conv)
                if True:
                    # Implicit convection
                    flux = calc_flux(u[d], u_conv)
                    eq -= dot(rho*u[d]*u_conv, grad(v[d]))*dx
                    eq += flux*jump(v[d])*dS
                else:
                    # Explicit convection
                    for fac, uc in ((2, sim.data['up']), (-1, sim.data['upp'])):
                        flux = calc_flux(uc[d], uc)
                        eq -= fac*dot(rho*uc[d]*uc, grad(v[d]))*dx
                        eq += fac*flux*jump(v[d])*dS
                
                # Calculate SIPG penalty 
                mpm = sim.multi_phase_model
                mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
                P = sim.data['Vu'].ufl_element().degree()
                penalty1 = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=2.0)
                penalty2 = define_penalty(mesh, P, mu_min, mu_max, boost_factor=30, exponent=2.0)
                #penalty = 3*mu_max*(P+1)*(P + ndim)/mesh.hmin()
                sim.log.info('\nDG SIPG penalties u%d:  %.4e  %.4e' % (d, penalty1, penalty2))
                penalty_dS = dolfin.Constant(penalty1)
                penalty_ds = dolfin.Constant(penalty2)
                
                # Diffusion:
                # -∇⋅μ[(∇u) + (∇u)^T]
                eq += mu*dot(grad(u[d]), grad(v[d]))*dx
                eq += sd*mu*dot(u.dx(d), grad(v[d]))*dx
                
                # Symmetric Interior Penalty method for -∇⋅μ∇u
                eq -= avg(dot(n, mu*grad(u[d])))*jump(v[d])*dS
                eq -= avg(dot(n, mu*grad(v[d])))*jump(u[d])*dS
                
                # Symmetric Interior Penalty method for -∇⋅μ(∇u)^T
                eq -= sd*avg(dot(n, mu*u.dx(d)))*jump(v[d])*dS
                eq -= sd*avg(dot(n, mu*v.dx(d)))*jump(u[d])*dS
                
                # Symmetric Interior Penalty coercivity term
                eq += penalty_dS*jump(u[d])*jump(v[d])*dS
                
                # Pressure
                # ∇p
                eq -= p*v[d].dx(d)*dx
                eq += avg(p)*jump(v[d])*n[d]('+')*dS
                
                # Body force (gravity)
                # ρ g
                eq -= rho*g[d]*v[d]*dx
                
                # Dirichlet boundary
                dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
                for dbc in dirichlet_bcs:
                    # Divergence free criterion
                    eq += q*(u[d] - dbc.func())*n[d]*dbc.ds()
                    
                    # SIPG for -∇⋅μ∇u 
                    eq -= mu*dot(n, grad(u[d]))*v[d]*dbc.ds()
                    eq -= mu*dot(n, grad(v[d]))*u[d]*dbc.ds()
                    eq += mu*dot(n, grad(v[d]))*dbc.func()*dbc.ds()
                    
                    # SIPG for -∇⋅μ(∇u)^T
                    eq -= sd*mu*dot(n, u.dx(d))*v[d]*dbc.ds()
                    eq -= sd*mu*dot(n, v.dx(d))*u[d]*dbc.ds()
                    eq += sd*mu*dot(n, v.dx(d))*dbc.func()*dbc.ds()
                    
                    # Weak Dirichlet
                    parts = (penalty_ds + abs(rho*dot(n, u_conv)))
                    eq += parts*(u[d] - dbc.func())*v[d]*dbc.ds()
                    
                    # Pressure
                    eq += p*v[d]*n[d]*dbc.ds()
                
                if d == 1:
                    print 'Penalty:', penalty1, penalty2, 'Flux:', self.flux_type
        
        self.eq_coupled = dolfin.system(eq)
    
    def coupled_boundary_conditions(self):
        """
        Convert boundary conditions from segregated to coupled function spaces
        """
        coupled_dirichlet_bcs = []
        for i, name in enumerate(self.subspace_names):
            if self.vel_is_discontinuous and name.startswith('u'):
                # Use weak BCs if the velocity is DG
                continue
            
            V = self.subspaces[i]
            bcs = self.simulation.data['dirichlet_bcs'].get(name, [])
            for bc in bcs:
                bc_new = bc.copy_and_change_function_space(V)
                coupled_dirichlet_bcs.append(bc_new)
        
        return coupled_dirichlet_bcs
    
    @timeit
    def update_convection(self, t, dt):
        """
        Update terms used to linearise and discretise the convective term
        """
        ndim = self.simulation.ndim
        data = self.simulation.data
        
        # Update convective velocity field components
        for d in range(ndim):
            uic = data['u_conv%d' % d]
            uip =  data['up%d' % d]
            uipp = data['upp%d' % d]
            
            if self.is_first_timestep:
                uic.vector()[:] = uip.vector()[:]
            else:
                # Backwards difference formulation - standard linear extrapolation 
                uic.vector()[:] = 2.0*uip.vector()[:] - 1.0*uipp.vector()[:]
        
        self.is_first_timestep = False
    
    @timeit
    def solve_coupled(self):
        """
        Solve the coupled equations
        """
        # Assemble the equation system
        a, L = self.eq_coupled
        A, b = dolfin.assemble_system(a, L)
        
        # Apply strong boundary conditions
        for dbc in self.dirichlet_bcs:
            dbc.apply(A, b)
        
        # Solve the equation system
        self.coupled_solver.solve(A, self.coupled_func.vector(), b)
        
        # Copy the data from the coupled solution vector and to the individual components
        up, _ = self.coupled_func.split()
        u, p = up.split(True)
        
        funcs = list(u.split(True)) + [p]
        for name, func in zip(self.subspace_names, funcs):
            self.simulation.data[name].vector()[:] = func.vector()
    
    @timeit
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation        
        sim.hooks.simulation_started()
        t = sim.time
        it = sim.timestep
        
        # Check if there are non-zero values in the upp vectors
        maxabs = 0
        for d in range(sim.ndim):
            this_maxabs = abs(sim.data['upp%d' % d].vector().array()).max()
            maxabs = max(maxabs, this_maxabs)
        has_upp_start_values = maxabs > 0 
        
        # Previous-previous values are provided so we can start up with second order time stepping 
        if has_upp_start_values:
            sim.log.info('Initial values for upp are found and used')
            self.is_first_timestep = False
            if self.timestepping_method == BDF:
                self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
        
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            
            # Check if the simulation is done
            if t+dt > tmax + 1e-6:
                break
            
            # Advance one time step
            it += 1
            t += dt
            self.simulation.data['dt'].assign(dt)
            self.simulation.hooks.new_timestep(it, t, dt)
            
            # Extrapolate the convecting velocity to the new time step
            self.update_convection(t, dt)
            
            # Solve for the new time step
            self.solve_coupled()
            
            # Move u -> up, up -> upp and prepare for the next time step
            for d in range(self.simulation.ndim):
                u_new = self.simulation.data['u%d' % d]
                up = self.simulation.data['up%d' % d]
                upp = self.simulation.data['upp%d' % d]
                upp.assign(up)
                up.assign(u_new)
            
            # Change time coefficient to second order
            if self.timestepping_method == BDF:
                self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
