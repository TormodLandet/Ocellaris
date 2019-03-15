# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import (
    verify_key,
    timeit,
    linear_solver_from_input,
    create_vector_functions,
    shift_fields,
    velocity_change,
    matmul,
)
from . import Solver, register_solver, BDM
from .simple_equations import SimpleEquations
from ..solver_parts import (
    VelocityBDMProjection,
    setup_hydrostatic_pressure,
    SlopeLimiterVelocity,
    before_simulation,
    after_timestep,
)


# Solvers - default values, can be changed in the input file
SOLVER_U_OPTIONS = {
    'use_ksp': True,
    'petsc_ksp_type': 'gmres',
    'petsc_pc_type': 'asm',
    'petsc_ksp_initial_guess_nonzero': True,
    'petsc_ksp_view': 'DISABLED',
    'inner_iter_rtol': [1e-10] * 3,
    'inner_iter_atol': [1e-15] * 3,
    'inner_iter_max_it': [100] * 3,
}
SOLVER_P_OPTIONS = {
    'use_ksp': True,
    'petsc_ksp_type': 'gmres',
    'petsc_pc_type': 'hypre',
    'petsc_pc_hypre_type': 'boomeramg',
    'petsc_ksp_initial_guess_nonzero': True,
    'petsc_ksp_view': 'DISABLED',
    'inner_iter_rtol': [1e-10] * 3,
    'inner_iter_atol': [1e-15] * 3,
    'inner_iter_max_it': [100] * 3,
}

# Equations - default values, can be changed in the input file
USE_STRESS_DIVERGENCE = False
USE_LAGRANGE_MULTIPLICATOR = False
USE_GRAD_P_FORM = False
USE_GRAD_Q_FORM = True
HYDROSTATIC_PRESSURE_CALCULATION_EVERY_TIMESTEP = False
INCOMPRESSIBILITY_FLUX_TYPE = 'central'

# PIMPLE iteration control
MAX_INNER_ITER = 10
MAX_PRESSURE_ITER = 3
ALLOWABLE_ERROR_INNER = 1e-10
SKIP_MOM_PRED = False
DEBUG_RHS = False

# Relaxation for normal iterations and for the last
# iteration in each time step
ALPHA_U = 0.7
ALPHA_P = 0.2
ALPHA_U_LAST = 1.0
ALPHA_P_LAST = 1.0

# Approximation of Ã
NUM_ELEMENTS_IN_BLOCK = 1
LUMP_DIAGONAL = False
MASS_ONLY_A_TILDE = False


@register_solver('PIMPLE')
class SolverPIMPLE(Solver):
    description = 'The PIMPLE algorithm'

    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the PIMPLE algorithm as implemented
        in OpenFOAM and partially described in the PhD thesis of Jasak (1996;
        the PISO loop only)

        Starting with the coupled Navier-Stokes saddle point system:

            | A  B |   | u |   | d |
            |      | . |   | = |   |                                     (1)
            | C  0 |   | p |   | e |

        where e is not allways zero since we use weak BCs for the normal velocity.

        (1) Momentum prediction: guess a pressure p* and then solve for u*

            A u* = d - B p*

            Optionally apply relaxation to u* here

        (2) Reformulate the equation for u* into an approximate equation for u**

            Ã u** = (Ã - A) u* - B p** + d

            Require continuity of u** by setting C u** = e and inserting the
            expression for u** from the approximate momentum equation into
            the continuity equation (Ã needs to be easily invertible / block diag)

            CÃ⁻¹B p** = C u* - e - CÃ⁻¹A u* + CÃ⁻¹ d

            Solve for p**. We may optionally apply relaxation to p** here

        (3) Update the velocity based on the new pressure

            u** = u* - Ã⁻¹A u* - Ã⁻¹B p** + Ã⁻¹d

        (4) Go back to step 2 and compute a new pressure correction if requested
            This is the "PISO-loop" of PIMPLE, steps 2 and 3

        (5) Go back to step (1) and complete the whole computation again if
            required for convergence. If this is done then the algorithm is using
            PIMPLE mode and relaxation should be applied
         """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()
        self.hydrostatic_pressure = setup_hydrostatic_pressure(simulation, needs_initial_value=True)
        ph_every_timestep = 'p_hydrostatic' in sim.data

        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])

        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)

        # Define weak forms
        matrices = SimpleEquations(
            simulation,
            use_stress_divergence_form=self.use_stress_divergence_form,
            use_grad_p_form=self.use_grad_p_form,
            use_grad_q_form=self.use_grad_q_form,
            use_lagrange_multiplicator=self.use_lagrange_multiplicator,
            include_hydrostatic_pressure=ph_every_timestep,
            incompressibility_flux_type=self.incompressibility_flux_type,
            num_elements_in_block=self.num_elements_in_block,
            lump_diagonal=self.lump_diagonal,
            a_tilde_is_mass=self.mass_only_A_tilde,
        )
        self.matrices = matrices

        # Slope limiter for the momentum equation velocity components
        self.slope_limiter = SlopeLimiterVelocity(sim, sim.data['u'], 'u', vel_w=sim.data['u_conv'])
        self.using_limiter = self.slope_limiter.active

        # Projection for the velocity
        self.velocity_postprocessor = None
        if self.velocity_postprocessing == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(
                sim, sim.data['u'], incompressibility_flux_type=self.incompressibility_flux_type
            )

        # Matrix and vector storage
        self.A = self.B = self.C = self.At = self.Atinv = self.D = self.E = None
        self.AtinvB = self.CAtinvB = self.AtinvA = self.CAtinvA = None

        # Store number of iterations
        self.niters_u = None
        self.niters_p = None

    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation

        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(
            self.simulation, 'solver/u', default_parameters=SOLVER_U_OPTIONS
        )
        self.pressure_solver = linear_solver_from_input(
            self.simulation, 'solver/p', default_parameters=SOLVER_P_OPTIONS
        )

        # Lagrange multiplicator or remove null space via PETSc
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value(
            'solver/use_lagrange_multiplicator', USE_LAGRANGE_MULTIPLICATOR, 'bool'
        )
        if self.use_lagrange_multiplicator:
            self.remove_null_space = False

        # No need for special treatment if the pressure is coupled via outlet BCs
        if sim.data['outlet_bcs']:
            self.remove_null_space = False
            self.use_lagrange_multiplicator = False

        # Control the form of the governing equations
        self.use_stress_divergence_form = sim.input.get_value(
            'solver/use_stress_divergence_form', USE_STRESS_DIVERGENCE, 'bool'
        )
        self.use_grad_p_form = sim.input.get_value(
            'solver/use_grad_p_form', USE_GRAD_P_FORM, 'bool'
        )
        self.use_grad_q_form = sim.input.get_value(
            'solver/use_grad_q_form', USE_GRAD_Q_FORM, 'bool'
        )
        self.incompressibility_flux_type = sim.input.get_value(
            'solver/incompressibility_flux_type', INCOMPRESSIBILITY_FLUX_TYPE, 'string'
        )
        assert sim.data['Vu'].ufl_element().family() == 'Discontinuous Lagrange'

        # Velocity post_processing
        self.velocity_postprocessing = sim.input.get_value(
            'solver/velocity_postprocessing', BDM, 'string'
        )
        verify_key(
            'velocity post processing', self.velocity_postprocessing, ('none', BDM), 'SIMPLE solver'
        )

        # Quasi-steady simulation input
        self.steady_velocity_eps = sim.input.get_value(
            'solver/steady_velocity_stopping_criterion', None, 'float'
        )
        self.is_steady = self.steady_velocity_eps is not None

        # How to approximate A_tilde
        self.num_elements_in_block = sim.input.get_value(
            'solver/num_elements_in_A_tilde_block', NUM_ELEMENTS_IN_BLOCK, 'int'
        )
        self.lump_diagonal = sim.input.get_value(
            'solver/lump_A_tilde_diagonal', LUMP_DIAGONAL, 'bool'
        )
        self.mass_only_A_tilde = sim.input.get_value(
            'solver/A_tilde_is_mass_only', MASS_ONLY_A_TILDE, 'bool'
        )

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
        sim.data['uvw_start'] = dolfin.Function(Vcoupled)
        sim.ndofs += Vcoupled.dim() + Vp.dim()

        # Create assigner to extract split function from uvw and vice versa
        self.assigner_split = dolfin.FunctionAssigner([Vu] * sim.ndim, Vcoupled)
        self.assigner_merge = dolfin.FunctionAssigner(Vcoupled, [Vu] * sim.ndim)

        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)

    @timeit
    def momentum_prediction(self):
        """
        Solve the momentum prediction equation
        """
        sim = self.simulation
        u_star = sim.data['uvw_star']
        p_star = sim.data['p']

        # Assemble the LHS matrices only the first inner iteration
        if self.inner_iteration == 1:
            A, A_tilde, A_tilde_inv, B, C = self.matrices.assemble_matrices()
            self.A = dolfin.as_backend_type(A)
            self.A_tilde = dolfin.as_backend_type(A_tilde)
            self.A_tilde_inv = dolfin.as_backend_type(A_tilde_inv)
            self.B = dolfin.as_backend_type(B)
            self.C = dolfin.as_backend_type(C)
            self.D = dolfin.as_backend_type(self.matrices.assemble_D())

        # Save the u_star we have at the beginning of the iteration so that we
        # can later compute the change in u_star after the end of the vel update
        sim.data['uvw_start'].assign(u_star)

        # Optionally skip the momentum prediction entirely
        if sim.input.get_value('solver/skip_momentum_predictor', SKIP_MOM_PRED, 'bool'):
            self.niters_u = 0
            return 0

        # The equation system to solve
        lhs = self.A
        rhs = self.D - self.B * p_star.vector()

        # Add optional under relaxation
        if self.last_inner_iter:
            alpha = sim.input.get_value('solver/relaxation_u_last_iter', ALPHA_U_LAST, 'float')
        else:
            alpha = sim.input.get_value('solver/relaxation_u', ALPHA_U, 'float')
        if alpha != 1.0:
            # Implicit under relaxation
            rhs = alpha * rhs + (1 - alpha) * lhs * u_star.vector()
            #relax = (1 - alpha) / alpha
            #lhs = dolfin.as_backend_type(lhs.copy())
            #lhs.axpy(relax, self.A_tilde, False)
            #lhs.apply('insert')
            #rhs.axpy(relax, self.A_tilde * u_star.vector())
            #rhs.apply('insert')

        # Solve the linearised convection-diffusion system
        self.niters_u = self.velocity_solver.inner_solve(
            lhs, u_star.vector(), rhs, in_iter=self.inner_iteration, co_iter=self.co_inner_iter
        )

    def piso_loop(self):
        """
        Perform pressure correction and velocity update iteratively.
        This is the "PISO" part of the PIMPLE algorithm
        """
        N = self.simulation.input.get_value('solver/num_pressure_corr', MAX_PRESSURE_ITER, 'int')
        for i in range(N):
            reassemble = self.inner_iteration == 1 and i == 0
            last_piso_iter = i == N - 1
            err_p, div_err = self.pressure_correction(reassemble, last_piso_iter)
            err_u = self.velocity_update()
            if N != 1 and not "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX":
                txt = (
                    '    Pressure correction %2d - err u* %10.3e - '
                    'err p %10.3e - div inp %10.3e - Kry. iters %3d'
                    % (i + 1, err_u, err_p, div_err, self.niters_p)
                )
                self.simulation.log.info(txt)

        # The change in u_star from before the momentum prediction and up to now
        u_star_old = self.simulation.data['uvw_start']
        u_star_new = self.simulation.data['uvw_star']
        err_u = (u_star_new.vector() - u_star_old.vector()).norm('l2')

        return err_p, err_u

    @timeit
    def pressure_correction(self, reassemble, last_piso_iter):
        """
        PIMPLE pressure correction
        """
        sim = self.simulation
        u_star = sim.data['uvw_star']
        p_star = sim.data['p']
        minus_p_hat = self.simulation.data['p_hat']

        # Assemble only once per time step
        if reassemble:
            self.E = dolfin.as_backend_type(self.matrices.assemble_E())

            # Compute LHS
            self.AtinvB = matmul(self.A_tilde_inv, self.B, self.AtinvB)
            self.CAtinvB = matmul(self.C, self.AtinvB, self.CAtinvB)

            # Needed for RHS
            self.AtinvA = matmul(self.A_tilde_inv, self.A, self.AtinvA)
            self.CAtinvA = matmul(self.C, self.AtinvA, self.CAtinvA)

        # Compute the residual divergence
        U = u_star.vector()
        div = self.C * U - self.E
        div_err = div.norm('l2')

        # The equation system
        lhs = self.CAtinvB
        rhs = div - self.CAtinvA * U + self.C * (self.A_tilde_inv * self.D)

        if DEBUG_RHS:
            # Quantify RHS contributions
            c0 = (self.C * U).norm('l2')
            c1 = (self.E).norm('l2')
            c2 = (self.CAtinvA * U).norm('l2')
            c3 = (self.C * (self.A_tilde_inv * self.D)).norm('l2')
            cT = max([c0, c1, c2, c3])
            sim.log.info(
                '                              Pressure RHS contributions:'
                '  %5.2f  %5.2f  %.2e    %5.2f  %5.2f    %.2e'
                % (c0 / cT, c1 / cT, div_err / cT, c2 / cT, c3 / cT, cT)
            )

        # Inform PETSc about the pressure null space
        if self.remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p_star.vector())
                null_vec[:] = 1
                null_vec *= 1 / null_vec.norm("l2")

                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])

            # Make sure the null space is set on the matrix
            if self.inner_iteration == 1:
                lhs.set_nullspace(self.pressure_null_space)

            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(rhs)

        # Solve for the new pressure correction
        minus_p_hat.assign(p_star)
        self.niters_p = self.pressure_solver.inner_solve(
            lhs, p_star.vector(), rhs, in_iter=self.inner_iteration, co_iter=self.co_inner_iter
        )

        # Compute change from last iteration
        minus_p_hat.vector().axpy(-1.0, p_star.vector())
        minus_p_hat.vector().apply('insert')

        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p_star.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1) * dx2)
            pavg = dolfin.assemble(p_star * dx2) / vol
            p_star.vector()[:] -= pavg

        # Explicit relaxation
        if self.last_inner_iter and last_piso_iter:
            alpha = sim.input.get_value('solver/relaxation_p_last_iter', ALPHA_P_LAST, 'float')
        else:
            alpha = sim.input.get_value('solver/relaxation_p', ALPHA_P, 'float')
        if alpha != 1:
            p_star.vector().axpy(1 - alpha, minus_p_hat.vector())
            p_star.vector().apply('insert')

        return minus_p_hat.vector().norm('l2'), div_err

    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """
        p = self.simulation.data['p']
        uvw = self.simulation.data['uvw_star']
        
        # Explicit velocity update
        #Ati = self.A_tilde_inv
        #AtiA = self.AtinvA
        #AtiB = self.AtinvB
        #delta_u = Ati * self.D - AtiA * uvw.vector() - AtiB * p.vector()

        # Explicit velocity update
        Ati, A, B, D = self.A_tilde_inv, self.A, self.B, self.D
        delta_u = Ati * (D - A * uvw.vector() - B * p.vector())
        
        uvw.vector().axpy(1.0, delta_u)
        uvw.vector().apply('insert')

        print(
            uvw.vector()[0],
            uvw.vector()[10000],
            uvw.vector()[20000],
            uvw.vector()[30000],
            uvw.vector()[40000],
            delta_u.norm('l2'),
        )

        if DEBUG_RHS:
            # Quantify RHS contributions
            c0 = (self.AtinvA * uvw.vector()).norm('l2')
            c1 = (self.AtinvB * p.vector()).norm('l2')
            c2 = (self.A_tilde_inv * self.D).norm('l2')
            cT = max([c0, c1, c2])
            self.simulation.log.info(
                '                              VelUp contributions:'
                '  %5.2f  %5.2f  %5.2f    %.2e' % (c0 / cT, c1 / cT, c2 / cT, cT)
            )

        return delta_u.norm('l2')

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
        change = velocity_change(
            u1=self.simulation.data['u'],
            u2=self.simulation.data['u_unlim'],
            ui_tmp=self.simulation.data['ui_tmp'],
        )

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

        with dolfin.Timer('Ocellaris run IPCS-A solver'):
            while True:
                # Get input values, these can possibly change over time
                dt = sim.input.get_value('time/dt', required_type='float')
                tmax = sim.input.get_value('time/tmax', required_type='float')
                num_inner_iter = sim.input.get_value('solver/num_inner_iter', MAX_INNER_ITER, 'int')
                allowable_error_inner = sim.input.get_value(
                    'solver/allowable_error_inner', ALLOWABLE_ERROR_INNER, 'float'
                )

                # Check if the simulation is done
                if t + dt > tmax + 1e-6:
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

                # Run inner iterations, PIMPLE loop
                self.inner_iteration = 1
                self.last_inner_iter = False
                while True:
                    if self.inner_iteration == num_inner_iter:
                        self.last_inner_iter = True

                    self.co_inner_iter = num_inner_iter - self.inner_iteration
                    self.momentum_prediction()
                    err_p, err_u = self.piso_loop()
                    sim.log.info(
                        '  PIMPLE iteration %3d - err u* %10.3e - err p %10.3e'
                        ' - Num Krylov iters - u %3d - p %3d'
                        % (self.inner_iteration, err_u, err_p, self.niters_u, self.niters_p)
                    )
                    self.inner_iteration += 1

                    if self.last_inner_iter:
                        break
                    elif err_u < allowable_error_inner:
                        self.last_inner_iter = True
                exit()

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
