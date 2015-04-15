from __future__ import division
import dolfin
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error, timeit, linear_solver_from_input
from . import Solver, register_solver, BDF, CRANK_NICOLSON
from .ipcs_equations import MomentumPredictionEquation, PressureCorrectionEquation, VelocityUpdateEquation


# Default values, can be changed in the input file
SOLVER_U = 'gmres'
PRECONDITIONER_U = 'additive_schwarz'
SOLVER_P = 'minres'
PRECONDITIONER_P = 'hypre_amg'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'relative_tolerance': 1e-10,
                     'absolute_tolerance': 1e-15}

# Implemented timestepping methods
TIMESTEPPING_METHODS = (BDF, CRANK_NICOLSON)


@register_solver('IPCS')
class SolverIPCS(Solver):
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
        
        # Define the momentum prediction equations
        self.eqs_mom_pred = []
        for d in range(sim.ndim):
            beta = self.convection_schemes[d].blending_function
            eq = MomentumPredictionEquation(simulation, d, beta, self.timestepping_method)
            self.eqs_mom_pred.append(eq)
        
        # Define the pressure correction equation
        self.eq_pressure = PressureCorrectionEquation(simulation)
        
        # Define the velocity update equations
        self.eqs_vel_upd = []
        for d in range(sim.ndim):
            eq = VelocityUpdateEquation(simulation, d)
            self.eqs_vel_upd.append(eq)
        
        # Storage for preassembled matrices
        self.Au = [None]*sim.ndim
        self.Ap = None
        self.Au_upd = None
        self.pressure_null_space = None
        
        # Store number of iterations
        self.niters_u = [None] * sim.ndim
        self.niters_p = None
        self.niters_u_upd = [None] * sim.ndim
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Test for PETSc
        if not dolfin.has_linear_algebra_backend("PETSc"):
            report_error('Missing PETSc',
                         'DOLFIN has not been configured with PETSc '
                         'which is needed by Ocellaris.')
        dolfin.parameters["linear_algebra_backend"] = "PETSc"
        
        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(self.simulation, 'solver/u', SOLVER_U,
                                                        PRECONDITIONER_U, None, KRYLOV_PARAMETERS)
        self.pressure_solver = linear_solver_from_input(self.simulation, 'solver/p', SOLVER_P,
                                                        PRECONDITIONER_P, None, KRYLOV_PARAMETERS)
        self.pressure_solver.parameters['preconditioner']['structure'] = 'same'
        self.u_upd_solver = linear_solver_from_input(self.simulation, 'solver/u_upd', SOLVER_U,
                                                     PRECONDITIONER_U, None, KRYLOV_PARAMETERS)
        
        # Get convection schemes for the velocity
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
    
    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation
        
        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        
        # Create velocity functions. Keep both component and vector forms
        u_list, up_list, upp_list, u_conv, u_star = [], [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['uppp%d' % d] = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            sim.data['u_star%d' % d] = us = dolfin.Function(Vu)
            u_list.append(u)
            up_list.append(up)
            upp_list.append(upp)
            u_conv.append(uc)
            u_star.append(us)
        sim.data['u'] = dolfin.as_vector(u_list)
        sim.data['up'] = dolfin.as_vector(up_list)
        sim.data['upp'] = dolfin.as_vector(upp_list)
        sim.data['u_conv'] = dolfin.as_vector(u_conv)
        sim.data['u_star'] = dolfin.as_vector(u_star)
        self.u_tmp = dolfin.Function(Vu)
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)
    
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
            elif self.timestepping_method == BDF:
                uic.vector()[:] = 2.0*uip.vector()[:] - 1.0*uipp.vector()[:]
            elif self.timestepping_method == CRANK_NICOLSON:
                # These two methods seem to give exactly the same results
                # Ingram (2013) claims that the first has better stability properties
                # in "A new linearly extrapolated Crank-Nicolson time-stepping scheme 
                # for the Navier-Stokes equations" 
                
                # Ingram's Crank-Nicolson extrapolation method
                #uippp = data['uppp%d' % d]
                #uic.vector()[:] = uip.vector()[:] + 0.5*uipp.vector()[:] - 0.5*uippp.vector()[:]
                
                # Standard Crank-Nicolson linear extrapolation method
                uic.vector()[:] = 1.5*uip.vector()[:] - 0.5*uipp.vector()[:]
        
        # Update the convection blending factors
        for cs in self.convection_schemes:
            cs.update(t, dt, data['u'])
        
        self.is_first_timestep = False
    
    @timeit
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        err = 0.0
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            self.u_tmp.vector()[:] = us.vector()
            
            dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            eq = self.eqs_mom_pred[d]
            
            if self.inner_iteration == 1:
                # Assemble the A matrix only the first inner iteration
                self.Au[d] = eq.assemble_lhs()
                self.velocity_solver.parameters['preconditioner']['structure'] = 'same_nonzero_pattern'
            else:
                self.velocity_solver.parameters['preconditioner']['structure'] = 'same'
            
            A = self.Au[d]
            b = eq.assemble_rhs()
            
            # Solve the advection equation
            family = us.element().family()
            if family == 'Lagrange':
                for dbc in dirichlet_bcs:
                    dbc.apply(A, b)
            
            self.niters_u[d] = self.velocity_solver.solve(A, us.vector(), b)
            
            self.u_tmp.vector().axpy(-1, us.vector())
            err += self.u_tmp.vector().norm('l2')
        return err
    
    @timeit
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
        """
        p = self.simulation.data['p']
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        
        # Assemble the A matrix only the first inner iteration
        if self.inner_iteration == 1:
            self.Ap = self.eq_pressure.assemble_lhs()
        
        # The equation system to solve
        A = self.Ap
        b = self.eq_pressure.assemble_rhs()
        
        # Apply strong boundary conditions
        family = p.element().family()
        if family == 'Lagrange':
            for dbc in dirichlet_bcs:
                dbc.apply(A, b)
        
        if not dirichlet_bcs:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p.vector())
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                
                # Create null space basis object and attach to Krylov solver
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])
                self.pressure_solver.set_nullspace(self.pressure_null_space)
            
            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(b)
        
        # Temporarily store the old pressure
        p_hat = self.simulation.data['p_hat']
        p_hat.vector().zero()
        p_hat.vector().axpy(-1, p.vector())
        
        # Solve for new pressure
        self.niters_p = self.pressure_solver.solve(A, p.vector(), b)
        
        # Calculate p_hat = p_new - p_old 
        p_hat.vector().axpy(1, p.vector())
        
        return p_hat.vector().norm('l2')
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """    
        for d in range(self.simulation.ndim):
            eq = self.eqs_vel_upd[d]
            
            if self.Au_upd is None:
                self.Au_upd = eq.assemble_lhs()
            
            A = self.Au_upd
            b = eq.assemble_rhs()
            u_new = self.simulation.data['u%d' % d]
            
            # Apply the Dirichlet boundary conditions
            family = u_new.element().family()
            if family == 'Lagrange':
                # Apply strong boundary conditions
                dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
                for bc in dirichlet_bcs:
                    bc.apply(A, b)
            
            self.niters_u_upd[d] = self.u_upd_solver.solve(A, u_new.vector(), b)
    
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
            num_inner_iter = sim.input.get_value('solver/num_inner_iter', 1, 'int')
            allowable_error_inner = sim.input.get_value('solver/allowable_error_inner', 0, 'float')
            
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
            
            # Run inner iterations
            for self.inner_iteration in xrange(1, num_inner_iter+1):
                err_u_star = self.momentum_prediction(t, dt)
                err_p = self.pressure_correction()
                
                # Information from solvers regarding number of iterations needed to solve linear system
                niters = ['%3d u%d' % (ni, d) for d, ni in enumerate(self.niters_u)]
                niters.append('%3d p' % self.niters_p)
                solver_info = ' - iters: %s' % ' '.join(niters)
                
                # Convergence estimates
                sim.log.info('  Inner iteration %3d - Linf* %10.3e - Linfp %10.3e%s'
                             % (self.inner_iteration, err_u_star, err_p, solver_info) + 
                             '  u0max %10.3e' % abs(sim.data['u_star0'].vector().array()).max())
                
                if err_u_star < allowable_error_inner:
                    break
            
            self.velocity_update()
            
            # Move u -> up, up -> upp and prepare for the next time step
            for d in range(self.simulation.ndim):
                u_new = self.simulation.data['u%d' % d]
                up = self.simulation.data['up%d' % d]
                upp = self.simulation.data['upp%d' % d]
                uppp = self.simulation.data['uppp%d' % d]
                uppp.assign(upp)
                upp.assign(up)
                up.assign(u_new)
            
            # Change time coefficient to second order
            if self.timestepping_method == BDF:
                self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
