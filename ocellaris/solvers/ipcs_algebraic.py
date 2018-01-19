import numpy.linalg
import dolfin
from ocellaris.utils import (verify_key, timeit, linear_solver_from_input,
                             create_vector_functions, shift_fields,
                             velocity_change, matmul, split_form_into_matrix)
from . import Solver, register_solver, BDM
from .coupled_equations import define_dg_equations
from ..solver_parts import (VelocityBDMProjection, setup_hydrostatic_pressure,
                            SlopeLimiterVelocity, before_simulation,
                            after_timestep)

# Solvers - default values, can be changed in the input file
SOLVER_U_OPTIONS = {'use_ksp': True,
                    'petsc_ksp_type': 'gmres', 
                    'petsc_pc_type': 'asm',
                    'petsc_ksp_initial_guess_nonzero': True,
                    'petsc_ksp_view': 'DISABLED',
                    'inner_iter_rtol': [1e-10] * 3,
                    'inner_iter_atol': [1e-15] * 3,
                    'inner_iter_max_it': [100] * 3}
SOLVER_P_OPTIONS = {'use_ksp': True,
                    'petsc_ksp_type': 'gmres',
                    'petsc_pc_type': 'hypre',
                    'petsc_pc_hypre_type': 'boomeramg',
                    'petsc_ksp_initial_guess_nonzero': True,
                    'petsc_ksp_view': 'DISABLED',
                    'inner_iter_rtol': [1e-10] * 3,
                    'inner_iter_atol': [1e-15] * 3,
                    'inner_iter_max_it': [100] * 3}
MAX_INNER_ITER = 10
ALLOWABLE_ERROR_INNER = 1e-10

# Equations - default values, can be changed in the input file
USE_STRESS_DIVERGENCE = False
USE_LAGRANGE_MULTIPLICATOR = False
USE_GRAD_P_FORM = False
USE_GRAD_Q_FORM = True
HYDROSTATIC_PRESSURE_CALCULATION_EVERY_TIMESTEP = False
INCOMPRESSIBILITY_FLUX_TYPE = 'central'


@register_solver('IPCS-A')
class SolverIPCSA(Solver):
    description = 'Incremental Pressure Correction Scheme (algebraic form)'
    
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on an algebraic Incremental Pressure correction scheme.
        
        Starting with the coupled Navier-Stokes saddle point system:
        
            | M+A  B |   | u |   | d |
            |        | . |   | = |   |                                     (1)
            | C    0 |   | p |   | e |
        
        where e is not allways zero since we use weak BCs for the normal velocity.
        
        (1) Momentum prediction: guess a pressure p* and then solve for u*
        
            (M + A) u* = d - B p*
        
        (2) Use the incompressibility constraint "C u = e" on "M(u - u*) = -B(p - p*)"
            (which can be added to (1) to get close to the correct momentum equation)
            the result is an equation for p (remember that M is block diagonal in DG)
            
            C M⁻¹ B p = C M⁻¹ B p* + Cu* - e
        
        (3) Repeat (1) and (2) until ||u-u*|| and/or ||p-p*|| are sufficiently close to
            zero, then update the velocity based on 
        
            u = u* - M⁻¹ B (p - p*)
        """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()
        self.hydrostatic_pressure = setup_hydrostatic_pressure(simulation, needs_initial_value=True)
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        
        # Define weak forms
        self.define_weak_forms()
        
        # Slope limiter for the momentum equation velocity components
        self.slope_limiter = SlopeLimiterVelocity(sim, sim.data['u'], 'u', vel_w=sim.data['u_conv'])
        self.using_limiter = self.slope_limiter.active
        
        # Projection for the velocity
        self.velocity_postprocessor = None
        if self.velocity_postprocessing == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(sim, sim.data['u'], 
                incompressibility_flux_type=self.incompressibility_flux_type)
        
        # Matrix and vector storage
        self.MplusA = self.B = self.C = self.M = self.Minv = self.D = self.E = None
        self.MinvB = self.CMinvB = None
        
        # Store number of iterations
        self.niters_u = None
        self.niters_p = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(self.simulation, 'solver/u',
                                                        default_parameters=SOLVER_U_OPTIONS)
        self.pressure_solver = linear_solver_from_input(self.simulation, 'solver/p',
                                                        default_parameters=SOLVER_P_OPTIONS)
        
        # Lagrange multiplicator or remove null space via PETSc
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value('solver/use_lagrange_multiplicator',
                                                              USE_LAGRANGE_MULTIPLICATOR, 'bool')
        if self.use_lagrange_multiplicator:
            self.remove_null_space = False
        
        # No need for special treatment if the pressure is coupled via outlet BCs
        if sim.data['outlet_bcs']:
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
        verify_key('velocity post processing', self.velocity_postprocessing, ('none', BDM), 'SIMPLE solver')
        
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
    
    def define_weak_forms(self):
        sim = self.simulation
        self.Vuvw = sim.data['uvw_star'].function_space()
        Vp = sim.data['Vp']
        
        # The trial and test functions in a coupled space (to be split)
        func_spaces = [self.Vuvw, Vp]
        e_mixed = dolfin.MixedElement([fs.ufl_element() for fs in func_spaces])
        Vcoupled = dolfin.FunctionSpace(sim.data['mesh'], e_mixed)
        tests = dolfin.TestFunctions(Vcoupled)
        trials = dolfin.TrialFunctions(Vcoupled)
        
        # Split into components
        v = dolfin.as_vector(tests[0][:])
        u = dolfin.as_vector(trials[0][:])
        q = tests[-1]
        p = trials[-1]
        lm_trial = lm_test = None
        
        # Define the full coupled form and split it into subforms depending
        # on the test and trial functions
        eq = define_dg_equations(u, v, p, q, lm_trial, lm_test, self.simulation,
                                 include_hydrostatic_pressure=self.hydrostatic_pressure.every_timestep,
                                 incompressibility_flux_type=self.incompressibility_flux_type,
                                 use_grad_q_form=self.use_grad_q_form,
                                 use_grad_p_form=self.use_grad_p_form,
                                 use_stress_divergence_form=self.use_stress_divergence_form)
        mat, vec = split_form_into_matrix(eq, Vcoupled, Vcoupled, check_zeros=True)
        
        # Check matrix and vector shapes and that the matrix is a saddle point matrix
        assert mat.shape == (2, 2)
        assert vec.shape == (2,)
        assert mat[-1,-1] is None, 'Found p-q coupling, this is not a saddle point system!'
        
        # Compile and store the forms
        self.eqA = dolfin.Form(mat[0, 0])
        self.eqB = dolfin.Form(mat[0, 1])
        self.eqC = dolfin.Form(mat[1, 0])
        self.eqD = dolfin.Form(vec[0])
        self.eqE = dolfin.Form(vec[1]) if vec[1] is not None else None
        
        # The mass matrix. Consistent with the implementation in define_dg_equations
        rho = sim.multi_phase_model.get_density(0)
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        eqM = rho*c1/dt*dolfin.dot(u, v)*dolfin.dx
        matM, _vecM = split_form_into_matrix(eqM, Vcoupled, Vcoupled, check_zeros=True)
        self.eqM = dolfin.Form(matM[0, 0])
    
    @timeit
    def momentum_prediction(self):
        """
        Solve the momentum prediction equation
        """
        sim = self.simulation
        u_star = sim.data['uvw_star']
        u_temp = sim.data['uvw_temp']
        p_star = sim.data['p']
        
        # Assemble only once per time step
        if self.inner_iteration == 1:
            self.MplusA = dolfin.as_backend_type(dolfin.assemble(self.eqA, tensor=self.MplusA))
            self.B = dolfin.as_backend_type(dolfin.assemble(self.eqB, tensor=self.B))
            self.D = dolfin.as_backend_type(dolfin.assemble(self.eqD, tensor=self.D))
        
        lhs = self.MplusA
        rhs = self.D - self.B * p_star.vector()
        
        self.project_rhs = False
        if self.project_rhs and self.velocity_postprocessor:
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
        
        # Solve the linearised convection-diffusion system
        u_temp.assign(u_star)
        self.niters_u = self.velocity_solver.inner_solve(lhs, u_star.vector(), rhs,
                                                         in_iter=self.inner_iteration,
                                                         co_iter=self.co_inner_iter)
        
        # Compute change from last iteration
        u_temp.vector().axpy(-1, u_star.vector())
        u_temp.vector().apply('insert')
        self._last_u_err = u_temp.vector().norm('l2')
        return self._last_u_err

    @timeit
    def compute_M_inverse(self):
        """
        Compute the inverse of the block diagonal mass matrix
        """
        M = self.M = dolfin.assemble(self.eqM, tensor=self.M)
        mesh = self.simulation.data['mesh']
        dm = self.Vuvw.dofmap()
        N = dm.cell_dofs(0).shape[0]
        Alocal = numpy.zeros((N, N), float)
        
        if self.Minv is None:
            self.Minv = M.copy()
            
        # Use mass matrix or block diagonal part of A
        self.splitting_approximation = 'mass'
        if self.splitting_approximation == 'mass':
            approx_A = M
        elif self.splitting_approximation == 'block_diag':
            approx_A = self.MplusA
        
        # Loop over cells and get the block diagonal parts (should be moved to C++)
        istart = M.local_range(0)[0]
        for cell in dolfin.cells(mesh, 'regular'):
            # Get global dofs
            dofs = dm.cell_dofs(cell.index()) + istart
            
            # Get block diagonal part of approx_A, invert it and insert into M⁻¹
            approx_A.get(Alocal, dofs, dofs)
            Alocal_inv = numpy.linalg.inv(Alocal)
            self.Minv.set(Alocal_inv, dofs, dofs)
        self.Minv.apply('insert')
        return self.Minv
    
    @timeit
    def pressure_correction(self):
        """
        Solve the Navier-Stokes equations on SIMPLE form
        (Semi-Implicit Method for Pressure-Linked Equations)
        """
        sim = self.simulation
        u_star = sim.data['uvw_star']
        p_star = sim.data['p']
        p_hat = self.simulation.data['p_hat']
        
        # Assemble only once per time step
        if self.inner_iteration == 1:
            self.C = dolfin.as_backend_type(dolfin.assemble(self.eqC, tensor=self.C))
            self.Minv = dolfin.as_backend_type(self.compute_M_inverse())
            if self.eqE is not None:
                self.E = dolfin.as_backend_type(dolfin.assemble(self.eqE, tensor=self.E))
            
            # Compute LHS
            self.MinvB = matmul(self.Minv, self.B, self.MinvB)
            self.CMinvB = matmul(self.C, self.MinvB, self.CMinvB)
        
        # The equation system
        lhs = self.CMinvB
        rhs = self.CMinvB * p_star.vector()
        rhs.axpy(1, self.C * u_star.vector())
        if self.eqE is not None:
            rhs.axpy(-1, self.E)
        rhs.apply('insert')
        
        # Inform PETSc about the pressure null space
        if self.remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p_star.vector())
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                
                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])
            
            # Make sure the null space is set on the matrix
            if self.inner_iteration == 1:
                lhs.set_nullspace(self.pressure_null_space)
            
            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(rhs)
        
        # Temporarily store the old pressure
        p_hat.vector().zero()
        p_hat.vector().axpy(-1, p_star.vector())
        
        # Solve for the new pressure correction
        self.niters_p = self.pressure_solver.inner_solve(lhs, p_star.vector(), rhs,
                                                         in_iter=self.inner_iteration,
                                                         co_iter=self.co_inner_iter)
        
        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p_star.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1)*dx2)
            pavg = dolfin.assemble(p_star*dx2)/vol
            p_star.vector()[:] -= pavg
        
        # Calculate p_hat = p_new - p_old 
        p_hat.vector().axpy(1, p_star.vector())
        
        return p_hat.vector().norm('l2')
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """
        p_hat = self.simulation.data['p_hat']
        uvw = self.simulation.data['uvw_star']
        uvw.vector().axpy(-1, self.MinvB * p_hat.vector())
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
    
        with dolfin.Timer('Ocellaris run IPCS-A solver'):    
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
                
                # Collect previous velocity components in coupled function
                self.assigner_merge.assign(sim.data['uvw_star'], list(sim.data['u']))
                
                # Run inner iterations
                self.inner_iteration = 1
                while self.inner_iteration <= num_inner_iter:
                    self.co_inner_iter = num_inner_iter - self.inner_iteration
                    err_u = self.momentum_prediction()
                    err_p = self.pressure_correction()
                    sim.log.info('  IPCS-A iteration %3d - err u* %10.3e - err p %10.3e'
                                 ' - Num Krylov iters - u %3d - p %3d' % (self.inner_iteration,
                                     err_u, err_p, self.niters_u, self.niters_p))
                    if err_u < allowable_error_inner:
                        break
                    
                    self.inner_iteration += 1
                    sim.flush()
                    
                    # Better initial guess for next solve
                    self.velocity_update()
                
                # Extract the separate velocity component functions
                self.assigner_split.assign(list(sim.data['u']), sim.data['uvw_star'])
                
                # Postprocess and limit velocity outside the inner iteration
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
