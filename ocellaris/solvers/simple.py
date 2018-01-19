import dolfin
from ocellaris.utils import (verify_key, timeit, linear_solver_from_input,
                             create_vector_functions, shift_fields,
                             velocity_change, matmul, dolfin_log_level)
from . import Solver, register_solver, BDM
from ..solver_parts import (VelocityBDMProjection, setup_hydrostatic_pressure,
                            SlopeLimiterVelocity, before_simulation,
                            after_timestep)
from .simple_equations import EQUATION_SUBTYPES


# Solvers - default values, can be changed in the input file
SOLVER_U = 'gmres'
PRECONDITIONER_U = 'additive_schwarz'
SOLVER_P = 'gmres'
PRECONDITIONER_P = 'hypre_amg'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'maximum_iterations': 100,
                     'relative_tolerance': 1e-10,
                     'absolute_tolerance': 1e-15,
                     'monitor_convergence': False,
                     'report': False}
MAX_INNER_ITER = 10
ALLOWABLE_ERROR_INNER = 1e-10

# Equations - default values, can be changed in the input file
EQUATION_SUBTYPE = 'Default'
USE_STRESS_DIVERGENCE = False
USE_LAGRANGE_MULTIPLICATOR = False
USE_GRAD_P_FORM = False
USE_GRAD_Q_FORM = True
HYDROSTATIC_PRESSURE_CALCULATION_EVERY_TIMESTEP = False
INCOMPRESSIBILITY_FLUX_TYPE = 'central'

ALPHA_U = 0.5
ALPHA_P = 0.5
LIMIT_INNER = False

NUM_ELEMENTS_IN_BLOCK = 0
LUMP_DIAGONAL = False
PROJECT_RHS = False


@register_solver('SIMPLE')
class SolverSIMPLE(Solver):
    def __init__(self, simulation, embedded=False):
        """
        A Navier-Stokes solver based on the Semi-Implicit Method for Pressure-Linked Equations (SIMPLE).
        
        Starting with the coupled Navier-Stokes saddle point system:
        
            | A  B |   | u |   | d |
            |      | . |   | = |   |                                     (1)
            | C  0 |   | p |   | 0 |
            
        Cannot solve for u since we do not know p, so we guess p* and get
        
            A u* = d - B p*
            C u* = e           <-- e is not necessarily zero since p* is not correct
            
        Subtracting from the real momentum equation and using
        
            u^ = u - u*
            p^ = p - p*
        
        we get
        
            A u^ = -B p^  
            C u^ = -e
        
        We can express u^ based on this
        
            u^ = - Ainv B p^                                             (8)
        
        and solve for p^ with Ã ≈ A (but easier to invert)
        
            C Ãinv B p^ = e                                              (9)
        
        We have to use under relaxation to update p and u (0 < α < 1)
        
            p = p* + α p^                                                (10)
        
        and for the velocity we use implicit under relaxation
        
            [(1-α)/α Ã + A] u* = d - B p* + (1-α)/α Ã u*_prev            (11)
        
        So
        
            1) Solve for u* using (11) with previous guesses u* and p*
            2) Find p^ using (9) and the divergence of u* from step 1
               to calculate e
            3) Update p using (10)
            4) Update u using (8)
            5) Check for convergence and possibly go back to 1 with new
               guesses for u* and p*
        
        Algorithm based on Klein, Kummer, Keil & Oberlack (2015) and DG discretsation
        based on Cockburn, Kanschat 
        """
        self.simulation = sim = simulation
        self.embedded = embedded
        self.read_input()
        if not embedded:
            self.create_functions()
            self.hydrostatic_pressure = setup_hydrostatic_pressure(simulation, needs_initial_value=True)
        ph_every_timestep = 'p_hydrostatic' in sim.data
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        
        # Get matrices
        Matrices = EQUATION_SUBTYPES[self.equation_subtype]
        matrices = Matrices(simulation,
                            use_stress_divergence_form=self.use_stress_divergence_form,
                            use_grad_p_form=self.use_grad_p_form,
                            use_grad_q_form=self.use_grad_q_form,
                            use_lagrange_multiplicator=self.use_lagrange_multiplicator,
                            include_hydrostatic_pressure=ph_every_timestep,
                            incompressibility_flux_type=self.incompressibility_flux_type,
                            num_elements_in_block=self.num_elements_in_block,
                            lump_diagonal=self.lump_diagonal)
        self.matrices = matrices
        
        # Slope limiter for the momentum equation velocity components
        self.slope_limiter = SlopeLimiterVelocity(sim, sim.data['u'], 'u', vel_w=sim.data['u_conv'])
        self.using_limiter = self.slope_limiter.active
        self.limit_inner_iterations = self.limit_inner_iterations and self.using_limiter
        
        # Projection for the velocity
        self.velocity_postprocessor = None
        if self.velocity_postprocessing == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(sim, sim.data['u'], 
                incompressibility_flux_type=self.incompressibility_flux_type)
        
        # Storage for preassembled matrices
        self.A = None
        self.A_tilde = None
        self.A_tilde_inv = None
        self.B = None
        self.C = None
        
        # Temporary matrices to store matrix matrix products
        self.mat_AinvB = None
        self.mat_CAinvB = None
        
        # Store number of iterations
        self.niters_u = None
        self.niters_p = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(self.simulation, 'solver/u', SOLVER_U,
                                                        PRECONDITIONER_U, None, KRYLOV_PARAMETERS)
        self.pressure_solver = linear_solver_from_input(self.simulation, 'solver/p', SOLVER_P,
                                                        PRECONDITIONER_P, None, KRYLOV_PARAMETERS)
        
        # Get under relaxation factors
        self.alpha_u = sim.input.get_value('solver/relaxation_u', ALPHA_U, 'float')
        self.alpha_p = sim.input.get_value('solver/relaxation_p', ALPHA_P, 'float')
        
        # Get the class to be used for the equation system assembly
        self.equation_subtype = sim.input.get_value('solver/equation_subtype', EQUATION_SUBTYPE, 'string')
        verify_key('equation sub-type', self.equation_subtype, EQUATION_SUBTYPES, 'SIMPLE solver')
        
        # Lagrange multiplicator or remove null space via PETSc
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value('solver/use_lagrange_multiplicator',
                                                              USE_LAGRANGE_MULTIPLICATOR, 'bool')
        if self.use_lagrange_multiplicator:
            self.remove_null_space = False
        
        # No need for special treatment if the pressure is set via Dirichlet conditions somewhere
        # No need for any tricks if the pressure is set via Dirichlet conditions somewhere
        if sim.data['dirichlet_bcs'].get('p', []) or sim.data['outlet_bcs']:
            self.remove_null_space = False
            self.use_lagrange_multiplicator = False
        
        # Control the form of the governing equations 
        self.use_stress_divergence_form = sim.input.get_value('solver/use_stress_divergence_form',
                                                              USE_STRESS_DIVERGENCE, 'bool')
        self.use_grad_p_form = sim.input.get_value('solver/use_grad_p_form', USE_GRAD_P_FORM, 'bool')
        self.use_grad_q_form = sim.input.get_value('solver/use_grad_q_form', USE_GRAD_Q_FORM, 'bool')
        self.incompressibility_flux_type = sim.input.get_value('solver/incompressibility_flux_type',
                                                               INCOMPRESSIBILITY_FLUX_TYPE, 'string')
        
        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Velocity post_processing
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing = sim.input.get_value('solver/velocity_postprocessing', default_postprocessing, 'string')
        self.project_rhs = sim.input.get_value('solver/project_rhs', PROJECT_RHS, 'bool')
        verify_key('velocity post processing', self.velocity_postprocessing, ('none', BDM), 'SIMPLE solver')
        
        # Quasi-steady simulation input
        self.steady_velocity_eps = sim.input.get_value('solver/steady_velocity_stopping_criterion',
                                                       None, 'float')
        self.is_steady = self.steady_velocity_eps is not None
        
        # Limiter inside inner iterations
        self.limit_inner_iterations = sim.input.get_value('solver/limit_inner_iterations', LIMIT_INNER, 'bool')
        
        # How to approximate A_tilde
        self.num_elements_in_block = sim.input.get_value('solver/num_elements_in_A_tilde_block', NUM_ELEMENTS_IN_BLOCK, 'int')
        self.lump_diagonal = sim.input.get_value('solver/lump_A_tilde_diagonal', LUMP_DIAGONAL, 'bool')
    
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
        
        def mom_assemble():
            """
            Assemble the linear systems
            """
            p_star = sim.data['p']
            alpha = self.alpha_u
            assert 0 < alpha
            relax = (1 - alpha) / alpha
            
            # Get coupled guess velocities
            uvw_star = sim.data['uvw_star']
            
            # Assemble the LHS matrices only the first inner iteration
            if self.inner_iteration == 1:
                self.A, self.A_tilde, self.A_tilde_inv, self.B, self.C = \
                    self.matrices.assemble_matrices()
            
            A = self.A
            A_tilde = self.A_tilde
            B = self.B
            D = self.matrices.assemble_D()
            
            if True:
                lhs = dolfin.as_backend_type(A.copy())
                lhs.axpy(relax, A_tilde, False)
                lhs.apply('insert')
                rhs = D
                rhs.axpy(-1.0, B * p_star.vector())
                rhs.axpy(relax, A_tilde * uvw_star.vector())
                rhs.apply('insert')
            
            else:
                lhs = dolfin.as_backend_type(A.copy())
                lhs._scale(1/alpha)
                lhs.apply('insert')
                rhs = D
                rhs.axpy(-1.0, B * p_star.vector())
                rhs.axpy(relax, A * uvw_star.vector())
                rhs.apply('insert')
            
            return lhs, rhs
        
        def mom_proj_rhs(rhs):
            """
            Project the RHS into BDM (embedded in DG)
            """
            if not self.project_rhs or not self.velocity_postprocessor:
                return
            print('PROJ PROJ PROJ PROJ PROJ PROJ PROJ PROJ PROJ!')
            
            if not hasattr(self, 'rhs_tmp'):
                # Setup RHS projection
                Vu = sim.data['Vu']
                funcs = [dolfin.Function(Vu) for _ in range(sim.ndim)]
                self.rhs_tmp = dolfin.as_vector(funcs)
                self.rhs_postprocessor = VelocityBDMProjection(sim, self.rhs_tmp,
                    incompressibility_flux_type=self.incompressibility_flux_type)
            
            self.assigner_split.assign(list(self.rhs_tmp), rhs)            
            self.rhs_postprocessor.run()
            self.assigner_merge.assign(rhs, list(self.rhs_tmp))
        
        def mom_solve(lhs, rhs):
            """
            Solve the linear systems
            """
            solver = self.velocity_solver
            if solver.is_iterative:
                if self.inner_iteration == 1:
                    solver.set_reuse_preconditioner(False)
                else:
                    solver.set_reuse_preconditioner(True)
                    #solver.parameters['maximum_iterations'] = 3
            
            u_star = sim.data['uvw_star']
            u_temp = sim.data['uvw_temp']
            u_temp.assign(u_star)
            
            with dolfin_log_level(dolfin.LogLevel.WARNING):
                self.niters_u = solver.solve(lhs, u_star.vector(), rhs)
            
            # Compute change from last iteration
            u_temp.vector().axpy(-1, u_star.vector())
            u_temp.vector().apply('insert')
            return u_temp.vector().norm('l2')
        
        # Assemble LHS and RHS, project RHS and finally solve the linear systems
        lhs, rhs = mom_assemble()
        mom_proj_rhs(rhs)
        err = mom_solve(lhs, rhs)
        
        return err
    
    @timeit
    def pressure_correction(self):
        """
        Solve the Navier-Stokes equations on SIMPLE form
        (Semi-Implicit Method for Pressure-Linked Equations)
        """
        solver = self.pressure_solver
        p_hat = self.simulation.data['p_hat']
        
        if solver.is_iterative:
            if self.inner_iteration == 1:
                solver.set_reuse_preconditioner(False)
            else:
                solver.set_reuse_preconditioner(True)
                #solver.parameters['maximum_iterations'] = 3
        
        # Compute the LHS = C⋅Ãinv⋅B
        if self.inner_iteration == 1 or self.limit_inner_iterations:
            C, Ainv, B = self.C, self.A_tilde_inv, self.B
            self.mat_AinvB = matmul(Ainv, B, self.mat_AinvB)
            self.mat_CAinvB = matmul(C, self.mat_AinvB, self.mat_CAinvB)
            self.LHS_pressure = dolfin.as_backend_type(self.mat_CAinvB.copy())
        LHS = self.LHS_pressure
        
        # Compute the divergence of u* and the rest of the right hand side
        uvw_star = self.simulation.data['uvw_star']
        RHS = self.matrices.assemble_E_star(uvw_star)
        
        # Inform PETSc about the null space
        if self.remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p_hat.vector())
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                
                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])
            
            # Make sure the null space is set on the matrix
            if self.inner_iteration == 1:
                LHS.set_nullspace(self.pressure_null_space)
            
            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(RHS)
        
        # Solve for the new pressure correction
        with dolfin_log_level(dolfin.LogLevel.WARNING):
            #print('STARTING PRESSURE KRYLOV SOLVER',
            #      'it=', self.simulation.timestep, 
            #      'iit=', self.inner_iteration)
            self.niters_p = solver.solve(LHS, p_hat.vector(), RHS)
        
        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here 
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p_hat.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1)*dx2)
            pavg = dolfin.assemble(p_hat*dx2)/vol
            p_hat.vector()[:] -= pavg
        
        # Calculate p = p* + α p^
        self.simulation.data['p'].vector().axpy(self.alpha_p, p_hat.vector())
        
        return p_hat.vector().norm('l2')
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """
        p_hat = self.simulation.data['p_hat']
        uvw = self.simulation.data['uvw_star']
        uvw.vector().axpy(-1, self.mat_AinvB * p_hat.vector())
        uvw.vector().apply('insert')
    
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
        div_dS_f, div_dx_f = self.simulation.solution_properties.divergences()
        div_dS = div_dS_f.vector().max()
        div_dx = div_dx_f.vector().max()
        return div_dS + div_dx
    
    @timeit.named('run SIMPLE solver')
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation
        sim.hooks.simulation_started()
        
        # Setup timestepping and initial convecting velocity
        before_simulation(sim)
        
        # Time loop
        t = sim.time
        it = sim.timestep
        
        # Give reasonable starting guesses for the solvers
        shift_fields(sim, ['up%d', 'u%d']) # get the initial u star
        
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            num_inner_iter = sim.input.get_value('solver/num_inner_iter',
                                                 MAX_INNER_ITER, 'int')
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
                # Collect previous velocity components in coupled function
                self.assigner_merge.assign(sim.data['uvw_star'], list(sim.data['u']))
                
                # Velocity and pressure inner solves
                err_u = self.momentum_prediction(t, dt)
                err_p = self.pressure_correction()
                self.velocity_update()
                solver_info = ' - Num Krylov iters - u %3d - p %3d' % (self.niters_u, self.niters_p)
                
                # Extract the separate velocity component functions
                self.assigner_split.assign(list(sim.data['u']), sim.data['uvw_star'])
                
                # Set u_conv equal to u
                #shift_fields(sim, ['u%d', 'u_conv%d'])
                
                # Postprocess and limit velocity inside the inner iteration
                if self.limit_inner_iterations:
                    self.postprocess_velocity()
                    shift_fields(sim, ['u%d', 'u_conv%d'])
                    err_lim = self.slope_limit_velocities()
                    solver_info += ' - err lim %10.3e' % err_lim
                    
                    # Divergence error, only relevant after self.postprocess_velocity() has been run
                    err_div = self.calculate_divergence_error()
                    solver_info += ' - err div %10.3e' % err_div
                
                # Convergence estimates
                sim.log.info('  SIMPLE iteration %3d - err u* %10.3e - err p %10.3e%s'
                             % (self.inner_iteration, err_u, err_p, solver_info))
                
                if err_u < allowable_error_inner:
                    break
                
                self.inner_iteration += 1
            
            # Postprocess and limit velocity outside the inner iteration
            if not self.limit_inner_iterations:
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
