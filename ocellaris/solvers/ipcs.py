from __future__ import division
import dolfin
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error, timeit, velocity_error_norm, create_krylov_solver
from . import Solver, register_solver
from .ipcs_equations import define_advection_problem, define_poisson_problem, define_penalty

# Default values, can be changed in the input file
SOLVER_U = 'gmres'
PRECONDITIONER_U = 'additive_schwarz'
SOLVER_P = 'minres'
PRECONDITIONER_P = 'hypre_amg'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True}

@register_solver('IPCS')
class SolverIPCS(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the pressure-velocity splitting
        scheme IPCS (Incremental Pressure Correction Scheme)
        """
        self.simulation = sim = simulation
        mesh = sim.data['mesh']
        
        # Test for PETSc
        if not dolfin.has_linear_algebra_backend("PETSc"):
            report_error('Missing PETSc',
                         'DOLFIN has not been configured with PETSc '
                         'which is needed by Ocellaris.')
        dolfin.parameters["linear_algebra_backend"] = "PETSc"
        
        # Function spaces
        Vu = simulation.data['Vu']
        Vp = simulation.data['Vp']
        
        # Solver for the velocity prediction
        velocity_solver_name = sim.input.get_value('solver/u/solver', SOLVER_U, 'string')
        velocity_preconditioner = sim.input.get_value('solver/u/preconditioner', PRECONDITIONER_U, 'string')
        velocity_solver_parameters = sim.input.get_value('solver/u/parameters', {}, 'dict')
        self.velocity_solver = create_krylov_solver(velocity_solver_name, velocity_preconditioner,
                                                    [KRYLOV_PARAMETERS, velocity_solver_parameters])
        
        # Make a solver for the pressure correction
        pressure_solver_name = sim.input.get_value('solver/p/solver', SOLVER_P, 'string')
        pressure_preconditioner = sim.input.get_value('solver/p/preconditioner', PRECONDITIONER_P, 'string')
        pressure_solver_parameters = sim.input.get_value('solver/p/parameters', {}, 'dict')
        self.pressure_solver = create_krylov_solver(pressure_solver_name, pressure_preconditioner,
                                                    [KRYLOV_PARAMETERS, pressure_solver_parameters])
        
        # Create velocity functions. Keep both component and vector forms
        uvec, upvec, uppvec, u_conv, u_star = [], [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            sim.data['u_star%d' % d] = us = dolfin.Function(Vu)
            uvec.append(u)
            upvec.append(up)
            uppvec.append(upp)
            u_conv.append(uc)
            u_star.append(us)
        sim.data['u'] = dolfin.as_vector(uvec)
        sim.data['up'] = dolfin.as_vector(upvec)
        sim.data['upp'] = dolfin.as_vector(uppvec)
        sim.data['u_conv'] = u_conv = dolfin.as_vector(u_conv)
        sim.data['u_star'] = u_star = dolfin.as_vector(u_star)
        
        # Create pressure function
        sim.data['p'] = p = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)
        
        # Mesh parameters
        n = dolfin.FacetNormal(mesh)
        
        # Physical properties
        rho = sim.data['rho']
        nu = sim.data['nu']
        g = sim.data['g']
        self.dt = dolfin.Constant(1.0)
        
        # Get convection schemes for the velocity
        conv_schemes = []
        for d in range(sim.ndim):
            conv_scheme_name = sim.input.get_value('convection/u/convection_scheme', required_type='string')
            conv_scheme = get_convection_scheme(conv_scheme_name)(simulation, 'u%d' % d)
            conv_schemes.append(conv_scheme)
        self.convection_schemes = conv_schemes
        
        # Coefficients for u, up and upp 
        self.time_coeffs = dolfin.Constant([1, -1, 0]) # First time step
        self.relaxation_value = 1.0
        self.relaxation = dolfin.Constant(self.relaxation_value)
        
        # Calculate SIPG penalties for elliptic DG solver
        Pu = Vu.ufl_element().degree()
        k_u = k_u_max = k_u_min = nu
        penalty_u = define_penalty(mesh, Pu, k_u_max, k_u_min)
        Pp = Vu.ufl_element().degree()
        k_p = k_p_max = k_p_min = 1/rho
        penalty_p = define_penalty(mesh, Pp, k_p_max, k_p_min)
        sim.log.info('\nDG SIPG penalties:\n  u: %.4e\n  p: %.4e' % (penalty_u, penalty_p))
        
        # Define the momentum prediction equations
        penalty = dolfin.Constant(penalty_u)
        self.eqs_mom_pred = []
        for d in range(sim.ndim):
            trial = dolfin.TrialFunction(Vu)
            test = dolfin.TestFunction(Vu)
            beta = conv_schemes[d].blending_function
            f = -1/rho*p.dx(d) + g[d]
            dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            neumann_bcs = self.simulation.data['neumann_bcs'].get('u%d' % d, [])
            a1, L1 = define_advection_problem(trial, test, upvec[d], uppvec[d],
                                              u_conv, n, beta, self.time_coeffs, self.dt,
                                              dirichlet_bcs)
            a2, L2 = define_poisson_problem(trial, test, k_u, f, n, penalty, dirichlet_bcs, neumann_bcs)
            eq = a1+a2, L1+L2
            self.eqs_mom_pred.append(eq)
        
        # Define the pressure correction equation
        trial = dolfin.TrialFunction(Vp)
        test = dolfin.TestFunction(Vp)
        penalty = dolfin.Constant(penalty_p)
        f = -self.time_coeffs[0]/self.dt * dolfin.nabla_div(u_star)
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        neumann_bcs = self.simulation.data['neumann_bcs'].get('p', [])
        self.eq_pressure = define_poisson_problem(trial, test, k_p, f, n, penalty, dirichlet_bcs, neumann_bcs)
        
        # For error norms in the convergence estimates
        elem = uvec[0].element()
        self.Vu_highp = dolfin.FunctionSpace(mesh, elem.family(), elem.degree() + 3)
        
        # Storage for preassembled matrices
        self.Au = [None] * sim.ndim 
        self.Ap = None
        
        # Store number of iterations
        self.niters_u = [None] * sim.ndim
        self.niters_p = None
    
    @timeit
    def update_convection(self, t, dt):
        """
        Update terms used to linearise and discretise the convective term
        """
        ndim = self.simulation.ndim
        data = self.simulation.data
        timestep = self.simulation.timestep
        
        # Update convective velocity field components
        for d in range(ndim):
            uic = data['u_conv%d' % d]
            uip =  data['up%d' % d]
            uipp = data['upp%d' % d]
            
            if timestep < 3:
                uic.vector()[:] = uip.vector()[:]
            else:
                uic.vector()[:] = 2*uip.vector()[:] - uipp.vector()[:]
            
            # Apply the Dirichlet boundary conditions
            if False:
                dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
                for bc in dirichlet_bcs:
                    bc.apply(uic.vector())
        
        # Update the convection blending factors
        for cs in self.convection_schemes:
            cs.update(t, dt, data['u'])
    
    @timeit
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            a, L = self.eqs_mom_pred[d]
            
            # Assemble the A matrix only the first inner iteration
            if self.inner_iteration == 1:
                A = dolfin.assemble(a)
                self.Au[d] = A
            else:
                A = self.Au[d]
            b = dolfin.assemble(L)
            
            # Solve the advection equation
            family = us.element().family()
            if family == 'Lagrange':
                for dbc in dirichlet_bcs:
                    dbc.apply(A, b)
            
            self.niters_u[d] = self.velocity_solver.solve(A, us.vector(), b)
    
    @timeit
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
        """
        p_hat = self.simulation.data['p_hat']
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        a, L = self.eq_pressure
        
        # Assemble the A matrix only the first inner iteration
        if self.inner_iteration == 1:
            A = dolfin.assemble(a)
            self.Ap = A
        else:
            A = self.Ap
        b = dolfin.assemble(L)
        
        # Apply strong boundary conditions
        family = p_hat.element().family()
        if family == 'Lagrange':
            for dbc in dirichlet_bcs:
                dbc.apply(A, b)
        
        if not dirichlet_bcs:
            # Create vector that spans the null space
            null_vec = dolfin.Vector(p_hat.vector())
            null_vec[:] = 1
            null_vec *= 1/null_vec.norm("l2")
            
            # Create null space basis object and attach to Krylov solver
            null_space = dolfin.VectorSpaceBasis([null_vec])
            self.pressure_solver.set_nullspace(null_space)
            
            # Orthogonalize b with respect to the null space
            null_space.orthogonalize(b)
        
        self.niters_p = self.pressure_solver.solve(A, p_hat.vector(), b)
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity estimates with the updated pressure
        field from the pressure correction equation
        """
        rho = self.simulation.data['rho']
        p_hat = self.simulation.data['p_hat']
        Vu = self.simulation.data['Vu']
        c1 = self.time_coeffs[0]
        
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            u_new = self.simulation.data['u%d' % d]
            
            # Update the velocity
            f = us - self.relaxation*self.dt/(c1*rho) * p_hat.dx(d)
            un = dolfin.project(f, Vu)
            u_new.assign(un)
            
            # Apply the Dirichlet boundary conditions
            family = u_new.element().family()
            if family == 'Lagrange':
                # Apply strong boundary conditions
                dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
                for bc in dirichlet_bcs:
                    bc.apply(u_new.vector())
    
    @timeit
    def pressure_update(self):
        """
        Update the pressure at the end of an inner iteration
        """
        p_hat = self.simulation.data['p_hat']
        p = self.simulation.data['p']
        p.vector()[:] += self.relaxation_value*p_hat.vector()[:] 
    
    @timeit
    def velocity_update_final(self):
        """
        Update the velocities at the end of the time step
        """
        for d in range(self.simulation.ndim):
            u_new = self.simulation.data['u%d' % d]
            up = self.simulation.data['up%d' % d]
            upp = self.simulation.data['upp%d' % d]
            upp.vector()[:] = up.vector()[:]
            up.vector()[:] = u_new.vector()[:]
    
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation        
        sim.hooks.simulation_started()
        t = 0
        it = 0
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            num_inner_iter = sim.input.get_value('solver/num_inner_iter', 1, 'int')
            allowable_error_inner = sim.input.get_value('solver/allowable_error_inner', 1e100, 'float')
            
            # Check if the simulation is done
            if t+dt > tmax + 1e-8:
                break
            
            # Advance one time step
            it += 1
            t += dt
            self.simulation.hooks.new_timestep(it, t, dt)
            self.dt.assign(dt)
            
            # Extrapolate the convecting velocity to the new time step
            self.update_convection(t, dt)
            
            for self.inner_iteration in xrange(1, num_inner_iter+1):
                self.momentum_prediction(t, dt)
                self.pressure_correction()
                self.velocity_update()
                self.pressure_update()
                
                # Convergence estimates in L2 norm
                if sim.input.get_value('solver/calculate_L2_norms', False, 'bool'):
                    L2s = velocity_error_norm(sim.data['u'], sim.data['u_star'], self.Vu_highp, 'L2')
                    L2c = velocity_error_norm(sim.data['u'], sim.data['u_conv'], self.Vu_highp, 'L2')
                    L2_info = ' - L2* %10.3e - L2c %10.3e' % (L2s, L2c)
                else:
                    L2_info = ''
                    
                # Solver information
                niters = ['%3d u%d' % (ni, d) for d, ni in enumerate(self.niters_u)]
                niters.append('%3d p' % self.niters_p)
                solver_info = ' - iters: %s' % ' '.join(niters)
                
                # Convergence estimates in L_infinity norm
                Linfs = velocity_error_norm(sim.data['u'], sim.data['u_star'], self.Vu_highp, 'Linf')
                Linfc = velocity_error_norm(sim.data['u'], sim.data['u_conv'], self.Vu_highp, 'Linf') 
                sim.log.info('  Inner iteration %3d - Linf* %10.3e - Linfc %10.3e%s%s'
                             % (self.inner_iteration, Linfs, Linfc, L2_info, solver_info))
                
                if Linfs < allowable_error_inner:
                    break
            
            self.velocity_update_final()
            
            # Change time coefficient to second order
            self.time_coeffs.assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
