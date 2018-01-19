import dolfin
from ocellaris.utils import (verify_key, timeit, linear_solver_from_input,
                             create_vector_functions, shift_fields,
                             velocity_change)
from . import Solver, register_solver, get_solver, BDM
from ..solver_parts import (VelocityBDMProjection, setup_hydrostatic_pressure,
                            SlopeLimiterVelocity, before_simulation,
                            after_timestep)
from .ipcs_equations import EQUATION_SUBTYPES
from .simple import dolfin_log_level


# Solvers - default values, can be changed in the input file
SOLVER_U = 'gmres'
PRECONDITIONER_U = 'additive_schwarz'
SOLVER_P = 'gmres'
PRECONDITIONER_P = 'hypre_amg'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'relative_tolerance': 1e-10,
                     'absolute_tolerance': 1e-15}
MAX_ITER_MOMENTUM = 1000
MAX_INNER_ITER = 10
ALLOWABLE_ERROR_INNER = 1e-10
NUM_SIMPLE_INNER_ITER = 0

# Equations - default values, can be changed in the input file
EQUATION_SUBTYPE = 'Default'
USE_STRESS_DIVERGENCE = False
USE_LAGRANGE_MULTIPLICATOR = False
USE_GRAD_P_FORM = False
HYDROSTATIC_PRESSURE_CALCULATION_EVERY_TIMESTEP = False
INCOMPRESSIBILITY_FLUX_TYPE = 'central'


@register_solver('IPCS')
class SolverIPCS(Solver):
    description = 'Incremental Pressure Correction Scheme (differential form)'
    
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on a pressure-velocity splitting scheme, 
        IPCS (Incremental Pressure Correction Scheme) on differential form,
        i.e., the equation for the pressure is formed from an elliptic weak
        form with given boundary conditions and jump stabilization terms.
        """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()
        self.hydrostatic_pressure = setup_hydrostatic_pressure(simulation, needs_initial_value=True)
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        
        # Get equations
        MomentumPredictionEquation, PressureCorrectionEquation, \
            VelocityUpdateEquation = EQUATION_SUBTYPES[self.equation_subtype]
        
        # Define the momentum prediction equations
        self.eqs_mom_pred = []
        self.eqs_mom_pred = MomentumPredictionEquation(simulation,
                                use_stress_divergence_form=self.use_stress_divergence_form,
                                use_grad_p_form=self.use_grad_p_form,
                                include_hydrostatic_pressure=self.hydrostatic_pressure.every_timestep)
        
        # Define the pressure correction equation
        self.eq_pressure = PressureCorrectionEquation(simulation,
                                use_lagrange_multiplicator=self.use_lagrange_multiplicator,
                                incompressibility_flux_type=self.incompressibility_flux_type)
        
        # Define the velocity update equations
        self.eqs_vel_upd = []
        for d in range(sim.ndim):
            eq = VelocityUpdateEquation(simulation, d)
            self.eqs_vel_upd.append(eq)
        
        # Slope limiter for the momentum equation velocity components
        self.slope_limiter = SlopeLimiterVelocity(sim, sim.data['u'], 'u', vel_w=sim.data['u_conv'])
        self.using_limiter = self.slope_limiter.active
        
        # Projection for the velocity
        self.velocity_postprocessor = None
        if self.velocity_postprocessing == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(sim, sim.data['u'],
                incompressibility_flux_type=self.incompressibility_flux_type)
        
        # Storage for preassembled matrices
        self.Au = [None]*sim.ndim
        self.Au_upd = None
        
        # Store number of iterations
        self.niters_u = [None] * sim.ndim
        self.niters_p = None
        self.niters_u_upd = [None] * sim.ndim
        
        # Storage for convergence checks
        self._error_cache = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(self.simulation, 'solver/u', SOLVER_U,
                                                        PRECONDITIONER_U, None, KRYLOV_PARAMETERS)
        self.pressure_solver = linear_solver_from_input(self.simulation, 'solver/p', SOLVER_P,
                                                        PRECONDITIONER_P, None, KRYLOV_PARAMETERS)
        #self.pressure_solver.parameters['preconditioner']['structure'] = 'same'
        
        # Velocity update can be performed with local solver for DG velocities
        self.use_local_solver_for_update = sim.input.get_value('solver/u_upd_local',
                                                               self.vel_is_discontinuous, 'bool')
        if self.use_local_solver_for_update:
            self.u_upd_solver = None # Will be set when LHS is ready
        else:
            self.u_upd_solver = linear_solver_from_input(self.simulation, 'solver/u_upd', SOLVER_U,
                                                         PRECONDITIONER_U, None, KRYLOV_PARAMETERS)
        
        # Get the class to be used for the equation system assembly
        self.equation_subtype = sim.input.get_value('solver/equation_subtype', EQUATION_SUBTYPE, 'string')
        verify_key('equation sub-type', self.equation_subtype, EQUATION_SUBTYPES, 'ipcs solver')
        
        # Lagrange multiplicator or remove null space via PETSc
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value('solver/use_lagrange_multiplicator',
                                                              USE_LAGRANGE_MULTIPLICATOR, 'bool')
        has_dirichlet = self.simulation.data['dirichlet_bcs'].get('p', []) or sim.data['outlet_bcs']
        if self.use_lagrange_multiplicator or has_dirichlet:
            self.remove_null_space = False
        
        # No need for special treatment if the pressure is set via Dirichlet conditions somewhere
        if has_dirichlet:
            self.use_lagrange_multiplicator = False
            self.remove_null_space = False
        
        # Control the form of the governing equations
        self.use_stress_divergence_form = sim.input.get_value('solver/use_stress_divergence_form',
                                                              USE_STRESS_DIVERGENCE, 'bool')
        self.use_grad_p_form = sim.input.get_value('solver/use_grad_p_form', USE_GRAD_P_FORM, 'bool')
        self.incompressibility_flux_type = sim.input.get_value('solver/incompressibility_flux_type',
                                                               INCOMPRESSIBILITY_FLUX_TYPE, 'string')
        
        # Velocity post_processing
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing = sim.input.get_value('solver/velocity_postprocessing', default_postprocessing, 'string')
        verify_key('velocity post processing', self.velocity_postprocessing, (None, BDM), 'ipcs solver')
        
        # Quasi-steady simulation input
        self.steady_velocity_eps = sim.input.get_value('solver/steady_velocity_stopping_criterion',
                                                       None, 'float')
        self.is_steady = self.steady_velocity_eps is not None
    
    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation
        
        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        
        # Create velocity functions on component and vector form
        create_vector_functions(sim, 'u', 'u%d', Vu)
        create_vector_functions(sim, 'up', 'up%d', Vu)
        create_vector_functions(sim, 'upp', 'upp%d', Vu)
        create_vector_functions(sim, 'u_conv', 'u_conv%d', Vu)
        create_vector_functions(sim, 'up_conv', 'up_conv%d', Vu)
        create_vector_functions(sim, 'upp_conv', 'upp_conv%d', Vu)
        create_vector_functions(sim, 'u_unlim', 'u_unlim%d', Vu)
        sim.data['ui_tmp'] = dolfin.Function(Vu)
        
        # Create coupled vector function
        ue = Vu.ufl_element()
        e_mixed = dolfin.MixedElement([ue] * sim.ndim)
        Vcoupled = dolfin.FunctionSpace(Vu.mesh(), e_mixed)
        sim.data['uvw_star'] = dolfin.Function(Vcoupled)
        sim.data['uvw_temp'] = dolfin.Function(Vcoupled)
        sim.ndofs += Vcoupled.dim() + Vp.dim()
        
        # Create assigner to extract split function from uvw and vice versa
        self.assigner_split = dolfin.FunctionAssigner([Vu] * sim.ndim, Vcoupled)
        self.assigner_merge = dolfin.FunctionAssigner(Vcoupled, [Vu] * sim.ndim)
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)
    
    @timeit
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        sim = self.simulation
        solver = self.velocity_solver
        
        # Collect previous velocity components in coupled function
        uvw_star = sim.data['uvw_star']
        self.assigner_merge.assign(uvw_star, list(sim.data['u']))
        uvw_temp = sim.data['uvw_temp']
        uvw_temp.assign(uvw_star)
        
        eq = self.eqs_mom_pred
        
        if self.inner_iteration == 1:
            # Assemble the A matrix only the first inner iteration
            self.Au = eq.assemble_lhs()
        
        A = self.Au
        b = eq.assemble_rhs()
        #solver.parameters['maximum_iterations'] = MAX_ITER_MOMENTUM
        with dolfin_log_level(dolfin.LogLevel.WARNING):
            self.niters_u = solver.solve(A, uvw_star.vector(), b)
        
        # Extract the separate velocity component functions
        self.assigner_split.assign(list(sim.data['u']), uvw_star)
        
        # Compute change from last iteration
        uvw_temp.vector().axpy(-1, uvw_star.vector())
        uvw_temp.vector().apply('insert')
        return uvw_temp.vector().norm('l2')
    
    @timeit
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
        """
        p = self.simulation.data['p']
        p_hat = self.simulation.data['p_hat']
        
        # Temporarily store the old pressure
        p_hat.vector().zero()
        p_hat.vector().axpy(-1, p.vector())
        
        # Assemble the A matrix only the first inner iteration
        assemble_A = self.inner_iteration == 1
        
        # Solve the pressure equation    
        self.niters_p = self.eq_pressure.solve(self.pressure_solver, 
                                               assemble_A,
                                               self.remove_null_space)
    
        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here 
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1)*dx2)
            pavg = dolfin.assemble(p*dx2)/vol
            p.vector()[:] -= pavg
        
        # Calculate p_hat = p_new - p_old 
        p_hat.vector().axpy(1, p.vector())
        
        return p_hat.vector().norm('l2')
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """
        if self.use_local_solver_for_update:
            # Element-wise projection
            if self.u_upd_solver is None:
                self.u_upd_solver = dolfin.LocalSolver(self.eqs_vel_upd[0].form_lhs)
                self.u_upd_solver.factorize()
            
            Vu = self.simulation.data['Vu']
            for d in range(self.simulation.ndim):
                eq = self.eqs_vel_upd[d]
                b = eq.assemble_rhs()
                u_new = self.simulation.data['u%d' % d]
                self.u_upd_solver.solve_local(u_new.vector(), b, Vu.dofmap())
                self.niters_u_upd[d] = 0
        
        else:
            # Global projection
            for d in range(self.simulation.ndim):
                eq = self.eqs_vel_upd
                
                if self.Au_upd is None:
                    self.Au_upd = eq.assemble_lhs()
                
                A = self.Au_upd
                b = eq.assemble_rhs()
                u_new = self.simulation.data['u%d' % d]
                
                self.niters_u_upd[d] = self.u_upd_solver.solve(A, u_new.vector(), b)
    
    @timeit
    def postprocess_velocity(self):
        """
        Apply a post-processing operator to the given velocity field
        """
        if self.velocity_postprocessor:
            self.velocity_postprocessor.run()
    
    @timeit
    def slope_limit_velocities(self):
        """
        Run the slope limiter
        """
        if not self.using_limiter:
            return 0
        
        # Store unlimited velocities and then run limiter
        shift_fields(self.simulation, ['u%d', 'u_unlim%d'])
        self.slope_limiter.run()
        
        # Measure the change in the field after limiting (l2 norm)
        change = velocity_change(u1=self.simulation.data['u'],
                                 u2=self.simulation.data['u_unlim'],
                                 ui_tmp=self.simulation.data['ui_tmp'])
        
        return change
    
    @timeit
    def calculate_divergence_error(self):
        """
        Check the convergence towards zero divergence. This is just for user output
        """
        sim = self.simulation
        
        if self._error_cache is None:
            dot, grad, jump, avg = dolfin.dot, dolfin.grad, dolfin.jump, dolfin.avg
            dx, dS, ds = dolfin.dx, dolfin.dS, dolfin.ds
            
            u_star = sim.data['u']
            mesh = u_star[0].function_space().mesh()
            V = dolfin.FunctionSpace(mesh, 'DG', 1)
            n = dolfin.FacetNormal(mesh)
            u = dolfin.TrialFunction(V)
            v = dolfin.TestFunction(V)
            
            a = u*v*dx
            L = dot(avg(u_star), n('+'))*jump(v)*dS \
                + dot(u_star, n)*v*ds \
                - dot(u_star, grad(v))*dx
            
            local_solver = dolfin.LocalSolver(a, L)
            error_func = dolfin.Function(V)
            
            self._error_cache = (local_solver, error_func)
        
        local_solver, error_func = self._error_cache
        local_solver.solve_local_rhs(error_func)
        err_div = max(abs(error_func.vector().min()),
                      abs(error_func.vector().max()))
        
        # HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK 
        if self.inner_iteration == 20:
            print('HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK')
            from matplotlib import pyplot
            
            fig = pyplot.figure(figsize=(10, 8))
            a = dolfin.plot(error_func, cmap='viridis', backend='matplotlib')
            pyplot.colorbar(a)
            fig.savefig('test_%f.png' % sim.time)
        # HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
                
        return err_div
    
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation
        simple = None
        sim.hooks.simulation_started()
        
        # Setup timestepping and initial convecting velocity
        before_simulation(sim)
        
        # Time loop
        t = sim.time
        it = sim.timestep
        
        # Give reasonable starting guesses for the solvers
        shift_fields(sim, ['up%d', 'u%d'])
        
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            num_inner_iter = sim.input.get_value('solver/num_inner_iter',
                                                 MAX_INNER_ITER, 'int')
            num_simple_steps = sim.input.get_value('solver/num_simple_inner_iter',
                                                   NUM_SIMPLE_INNER_ITER, 'int')
            allowable_error_inner = sim.input.get_value('solver/allowable_error_inner',
                                                        ALLOWABLE_ERROR_INNER, 'float')
            
            # Check if the simulation is done
            if t+dt > tmax + 1e-6:
                break
            
            # Advance one time step
            it += 1
            t += dt
            self.simulation.data['dt'].assign(dt)
            self.simulation.hooks.new_timestep(it, t, dt)
            
            # Calculate the hydrostatic pressure when the density is not constant
            self.hydrostatic_pressure.update()
            
            # Run inner iterations
            self.inner_iteration = 1
            while self.inner_iteration <= num_inner_iter:
                err_u_star = self.momentum_prediction(t, dt)
                err_p = self.pressure_correction()
                
                # Information from solvers regarding number of iterations needed to solve linear system
                solver_info = ' - Num Krylov iters - u %3d - p %3d' % (self.niters_u, self.niters_p)
                
                # Get max u_star
                ustarmax = 0
                for d in range(sim.ndim):
                    thismax = abs(sim.data['u%d' % d].vector().get_local()).max()
                    ustarmax = max(thismax, ustarmax)
                ustarmax = dolfin.MPI.max(dolfin.MPI.comm_world, float(ustarmax))
                
                # Get the divergence error
                err_div = self.calculate_divergence_error()
                
                # Convergence estimates
                sim.log.info('  IPCS iteration %5d - err u* %10.3e - err p %10.3e%s  ui*max %10.3e'
                             % (self.inner_iteration, err_u_star, err_p, solver_info,  ustarmax)
                             + ' err div %10.3e' % err_div)
                
                if err_u_star < allowable_error_inner:
                    break
                
                self.inner_iteration += 1
            
            self.velocity_update()
            
            if num_simple_steps > 0:
                # Run SIMPLE iteration at end of the inner iterations to eliminate the divergence
                if simple is None:
                    simple = get_solver('SIMPLE')(sim, embedded=True)
                self.assigner_merge.assign(sim.data['uvw_star'], list(sim.data['u']))
                for simple_inner in range(num_simple_steps):
                    simple.inner_iteration = simple_inner + 1
                    err_u = simple.momentum_prediction(t, dt)
                    err_p = simple.pressure_correction()
                    simple.velocity_update()
                    solver_info = ' - Num Krylov iters - u %3d - p %3d' % (simple.niters_u, simple.niters_p)
                    sim.log.info('  SIMPLE iteration %3d - err u* %10.3e - err p %10.3e%s'
                                 % (simple.inner_iteration, err_u, err_p, solver_info))
                    if err_u_star < allowable_error_inner:
                        break
                self.assigner_split.assign(list(sim.data['u']), sim.data['uvw_star'])
            
            self.postprocess_velocity()
            shift_fields(sim, ['u%d', 'u_conv%d'])
            if self.using_limiter:
                self.slope_limit_velocities()
            
            # Move u -> up, up -> upp and prepare for the next time step
            vel_diff = after_timestep(sim, self.is_steady)
            
            # Stop steady state simulation if convergence has been reached
            if self.is_steady:
                vel_diff = dolfin.MPI.max(dolfin.MPI.comm_world, float(vel_diff))
                sim.reporting.report_timestep_value('max(ui_new-ui_prev)', vel_diff)                
                if vel_diff < self.steady_velocity_eps:
                    sim.log.info('Stopping simulation, steady state achieved')
                    sim.input.set_value('time/tmax', t)
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
