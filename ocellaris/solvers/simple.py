# encoding: utf8
from __future__ import division
import numpy
import dolfin
from ocellaris.utils import verify_key, timeit, linear_solver_from_input
from . import Solver, register_solver, BDM
from ..solver_parts import VelocityBDMProjection, HydrostaticPressure, SlopeLimiterVelocity
from .simple_equations import EQUATION_SUBTYPES, create_block_matrix


# Solvers - default values, can be changed in the input file
SOLVER_U = 'gmres'
PRECONDITIONER_U = 'additive_schwarz'
SOLVER_P = 'gmres'
PRECONDITIONER_P = 'hypre_amg'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'maximum_iterations': 10,
                     'relative_tolerance': 1e-10,
                     'absolute_tolerance': 1e-15}
MAX_INNER_ITER = 100
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


@register_solver('SIMPLE')
class SolverSIMPLE(Solver):
    def __init__(self, simulation):
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
        self.read_input()
        self.create_functions()
        self.setup_hydrostatic_pressure_calculations()
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        self.is_first_timestep = True
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        
        # Get matrices
        Matrices = EQUATION_SUBTYPES[self.equation_subtype]
        matrices = Matrices(simulation,
                            use_stress_divergence_form=self.use_stress_divergence_form,
                            use_grad_p_form=self.use_grad_p_form,
                            use_grad_q_form=self.use_grad_q_form,
                            use_lagrange_multiplicator=self.use_lagrange_multiplicator,
                            include_hydrostatic_pressure=self.hydrostatic_pressure_correction,
                            incompressibility_flux_type=self.incompressibility_flux_type)
        self.matrices = matrices
        
        # Slope limiter for the momentum equation velocity components
        self.slope_limiter = SlopeLimiterVelocity(sim, sim.data['u'], 'u')
        self.using_limiter = self.slope_limiter.active
        
        # Projection for the velocity
        self.velocity_postprocessor = None
        if self.velocity_postprocessing == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(sim, sim.data['u'], 
                incompressibility_flux_type=self.incompressibility_flux_type)
        
        # Storage for preassembled matrices
        self.A = [None]*sim.ndim
        self.A_tilde = [None]*sim.ndim
        self.A_tilde_inv = [None]*sim.ndim
        self.B = [None]*sim.ndim
        self.C = [None]*sim.ndim
        
        # Temporary matrices to store matrix matrix products
        self.mat_AinvB = [None]*sim.ndim
        self.mat_CAinvB = [None]*sim.ndim
        
        # Store number of iterations
        self.niters_u = [None] * sim.ndim
        self.niters_p = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        assert self.vel_is_discontinuous
        
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
        self.use_grad_q_form = sim.input.get_value('solver/use_grad_q_form', USE_GRAD_Q_FORM, 'bool')
        self.incompressibility_flux_type = sim.input.get_value('solver/incompressibility_flux_type',
                                                               INCOMPRESSIBILITY_FLUX_TYPE, 'string')
        
        # Velocity post_processing
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing = sim.input.get_value('solver/velocity_postprocessing', default_postprocessing, 'string')
        verify_key('velocity post processing', self.velocity_postprocessing, ('none', BDM), 'ipcs solver')
        
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
        
        # Create velocity functions. Keep both component and vector forms
        u_list, up_list, upp_list, u_conv = [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            sim.data['u_unlim%d' % d] = dolfin.Function(Vu)
            u_list.append(u)
            up_list.append(up)
            upp_list.append(upp)
            u_conv.append(uc)
        sim.data['u'] = dolfin.as_vector(u_list)
        sim.data['up'] = dolfin.as_vector(up_list)
        sim.data['upp'] = dolfin.as_vector(upp_list)
        sim.data['u_conv'] = dolfin.as_vector(u_conv)
        self.u_tmp = dolfin.Function(Vu)
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)
        
        # If gravity is nonzero we create a separate hydrostatic pressure field
        if any(gi != 0 for gi in sim.data['g'].values()):
            # Hydrostatic pressure is always CG
            Pp = Vp.ufl_element().degree()
            Vph = dolfin.FunctionSpace(sim.data['mesh'], 'CG', Pp)
            sim.data['p_hydrostatic'] = dolfin.Function(Vph)
    
    def setup_hydrostatic_pressure_calculations(self):
        """
        We can initialize the pressure field to be hydrostatic at the first time step,
        or we can calculate the hydrostatic pressure as its own pressure field every
        time step such that the PressureCorrectionEquation only solves for the dynamic
        pressure
        """
        sim = self.simulation
        
        # No need for hydrostatic pressure if g is zero
        g = sim.data['g']
        self.hydrostatic_pressure_correction = False
        if all(gi == 0 for gi in g.values()):
            return
        return
        
        # Get the input needed to calculate p_hydrostatic
        rho = sim.data['rho']
        sky_location = sim.input.get_value('multiphase_solver/sky_location', required_type='float')
        self.ph_every_timestep = sim.input.get_value('solver/hydrostatic_pressure_calculation_every_timestep', 
                                                     HYDROSTATIC_PRESSURE_CALCULATION_EVERY_TIMESTEP,
                                                     required_type='float')
        
        # Helper class to calculate the hydrostatic pressure distribution
        ph = sim.data['p_hydrostatic']
        self.hydrostatic_pressure = HydrostaticPressure(rho, g, ph, sky_location)
        self.hydrostatic_pressure_correction = True
    
    def update_hydrostatic_pressure(self):
        """
        Update the hydrostatic pressure field
        (when the density is not constant)
        """
        if not self.hydrostatic_pressure_correction:
            return
        
        self.hydrostatic_pressure.update()
            
        if not self.ph_every_timestep:
            # Correct the pressure only now, at the begining
            sim = self.simulation
            p = sim.data['p']
            if p.vector().max() == p.vector().min() == 0.0:
                sim.log.info('Initial pressure field is identically zero, initializing to hydrostatic')
                self.hydrostatic_pressure.update()
                p.assign(dolfin.interpolate(sim.data['p_hydrostatic'], p.function_space()))
            del sim.data['p_hydrostatic']
            self.hydrostatic_pressure_correction = False
    
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
                uic.assign(uip)
            else:
                uic.vector().zero()
                uic.vector().axpy(2.0, uip.vector())
                uic.vector().axpy(-1.0, uipp.vector())
        self.is_first_timestep = False
    
    @timeit
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        solver = self.velocity_solver
        p_star = self.simulation.data['p']
        
        if solver.is_iterative:
            if self.inner_iteration == 1:
                solver.set_reuse_preconditioner(False)
            else:
                solver.set_reuse_preconditioner(True)
                #solver.parameters['maximum_iterations'] = 3
        
        err = 0.0
        for d in range(self.simulation.ndim):
            u_star = self.simulation.data['u%d' % d]
            self.u_tmp.assign(u_star)
            
            # Assemble the LHS matrices only the first inner iteration
            if self.inner_iteration == 1:
                (self.A[d], self.A_tilde[d], self.A_tilde_inv[d],
                 self.B[d], self.C[d]) = self.matrices.assemble_matrices(d)
            
            A = self.A[d]
            A_tilde = self.A_tilde[d]
            B = self.B[d]
            D = self.matrices.assemble_D(d)
            alpha = self.alpha_u
            assert 0 < alpha
            relax = (1 - alpha) / alpha
            
            if True:
                LHS = dolfin.as_backend_type(A.copy())
                LHS.axpy(relax, A_tilde, False)
                LHS.apply('insert')
                RHS = D
                RHS.axpy(-1.0, B * p_star.vector())
                RHS.axpy(relax, A_tilde * u_star.vector())
                RHS.apply('insert')
            
            else:
                LHS = dolfin.as_backend_type(A.copy())
                LHS._scale(1/alpha)
                LHS.apply('insert')
                RHS = D
                RHS.axpy(-1.0, B * p_star.vector())
                RHS.axpy(relax, A * u_star.vector())
                RHS.apply('insert')
            
            with dolfin_log_level(dolfin.ERROR):
                self.niters_u[d] = solver.solve(LHS, u_star.vector(), RHS)
            
            self.u_tmp.vector().axpy(-1, u_star.vector())
            self.u_tmp.vector().apply('insert')
            err += self.u_tmp.vector().norm('l2')
        
        return err
    
    @timeit
    def slope_limit_velocities(self):
        """
        Run the slope limiter and assemble the limiter matrix
        """
        if not self.using_limiter:
            return 0
        
        for d in range(self.simulation.ndim):
            u_star = self.simulation.data['u%d' % d]
            u_unlim = self.simulation.data['u_unlim%d' % d]
            u_unlim.assign(u_star)
        
        self.slope_limiter.run()
        
        change = 0
        for d in range(self.simulation.ndim):
            u_star = self.simulation.data['u%d' % d]
            u_unlim = self.simulation.data['u_unlim%d' % d]
            self.u_tmp.assign(u_unlim)
            self.u_tmp.vector().axpy(-1, u_star.vector())
            self.u_tmp.vector().apply('insert')
            change += self.u_tmp.vector().norm('l2')
        
        return change
    
    @timeit
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
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
        if self.inner_iteration == 1 or self.using_limiter:
            LHS = 0
            for d in range(self.simulation.ndim):
                C, Ainv, B = self.C[d], self.A_tilde_inv[d], self.B[d]
                self.mat_AinvB[d] = matmul(Ainv, B, self.mat_AinvB[d])
                self.mat_CAinvB[d] = matmul(C, self.mat_AinvB[d], self.mat_CAinvB[d])
                
                if LHS == 0:
                    LHS = dolfin.as_backend_type(self.mat_CAinvB[d].copy())
                else:
                    LHS.axpy(1.0, self.mat_CAinvB[d], True)
                LHS.apply('insert')
            self.LHS_pressure = LHS
        else:
            LHS = self.LHS_pressure
        
        # Compute the divergence of u* and the rest of the right hand side
        u_star = self.simulation.data['u']
        RHS = self.matrices.assemble_E_star(u_star)
        
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
        with dolfin_log_level(dolfin.ERROR):
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
        
        for d in range(self.simulation.ndim):
            u = self.simulation.data['u%d' % d]
            u.vector().axpy(-1, self.mat_AinvB[d] * p_hat.vector())
    
    @timeit
    def postprocess_velocity(self):
        """
        Apply a post-processing operator to the given velocity field
        """
        if self.velocity_postprocessor:
            self.velocity_postprocessor.run()
            
    @timeit
    def calculate_divergence_error(self):
        """
        Check the convergence towards zero divergence. This is just for user output
        """
        div_dS_f, div_dx_f = self.simulation.solution_properties.divergences()
        div_dS = div_dS_f.vector().max()
        div_dx = div_dx_f.vector().max()
        return div_dS + div_dx
    
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
            this_maxabs = abs(sim.data['upp%d' % d].vector().get_local()).max()
            maxabs = max(maxabs, this_maxabs)
        maxabs = dolfin.MPI.max(dolfin.mpi_comm_world(), float(maxabs))
        has_upp_start_values = maxabs > 0
        
        # Previous-previous values are provided so we can start up with second order time stepping 
        if has_upp_start_values:
            sim.log.info('Initial values for upp are found and used')
            self.is_first_timestep = False
            self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
        
        # Give reasonable starting guesses for the solvers
        for d in range(sim.ndim):
            up = self.simulation.data['up%d' % d]
            u_new = self.simulation.data['u%d' % d]
            u_new.assign(up)
        
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
            
            # Extrapolate the convecting velocity to the new time step
            self.update_convection(t, dt)
            
            # Calculate the hydrostatic pressure when the density is not constant
            self.update_hydrostatic_pressure()
            
            # Run inner iterations
            self.inner_iteration = 1
            while self.inner_iteration <= num_inner_iter:
                err_u = self.momentum_prediction(t, dt)
                err_p = self.pressure_correction()
                
                self.velocity_update()
                self.postprocess_velocity()
                
                err_lim = self.slope_limit_velocities()
                
                # Information from solvers regarding number of iterations needed to solve linear system
                niters = ['%3d u%d' % (ni, d) for d, ni in enumerate(self.niters_u)]
                niters.append('%3d p' % self.niters_p)
                solver_info = ' - iters: %s' % ' '.join(niters)
                
                if self.using_limiter:
                    solver_info += ' - err lim %10.3e' % err_lim
                
                # Divergence error
                err_div = self.calculate_divergence_error()
                solver_info += ' - err div %10.3e' % err_div
                
                # Convergence estimates
                sim.log.info('  Inner iteration %3d - err u* %10.3e - err p %10.3e%s'
                             % (self.inner_iteration, err_u, err_p, solver_info))
                
                if err_u + err_p < allowable_error_inner:
                    break
                
                self.inner_iteration += 1
            
            #self.postprocess_velocity()
            
            # Move u -> up, up -> upp and prepare for the next time step
            vel_diff = 0
            for d in range(sim.ndim):
                u_new = sim.data['u%d' % d]
                up = sim.data['up%d' % d]
                upp = sim.data['upp%d' % d]
                
                if self.is_steady:
                    diff = abs(u_new.vector().get_local() - up.vector().get_local()).max() 
                    vel_diff = max(vel_diff, diff)
                
                upp.assign(up)
                up.assign(u_new)
            
            # Change time coefficient to second order
            sim.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Stop steady state simulation if convergence has been reached
            if self.is_steady:
                vel_diff = dolfin.MPI.max(dolfin.mpi_comm_world(), float(vel_diff))
                sim.reporting.report_timestep_value('max(ui_new-ui_prev)', vel_diff)                
                if vel_diff < self.steady_velocity_eps:
                    sim.log.info('Stopping simulation, steady state achieved')
                    sim.input.set_value('time/tmax', t)
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)


def matmul(A, B, out):
    """
    A B (and potentially out) must be PETScMatrix
    The matrix out must be the result of a prior matmul
    call with the same sparsity patterns in A and B
    """
    A = A.mat()
    B = B.mat()
    if out is not None:
        A.matMult(B, out.mat())
        C = out
    else:
        Cmat = A.matMult(B)
        C = dolfin.PETScMatrix(Cmat)
    
    C.apply('insert')
    return C


from contextlib import contextmanager
@contextmanager
def dolfin_log_level(level):
    old_level = dolfin.get_log_level()
    dolfin.set_log_level(level)
    yield
    dolfin.set_log_level(old_level)
