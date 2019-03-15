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
from ..solver_parts import (
    VelocityBDMProjection,
    setup_hydrostatic_pressure,
    SlopeLimiterVelocity,
    before_simulation,
    after_timestep,
)
from .simple_equations import EQUATION_SUBTYPES


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

SOLVER_SIMPLE = 'SIMPLE'
SOLVER_PISO = 'PISO'
ALPHA_U = 0.7
ALPHA_P = 0.9

NUM_ELEMENTS_IN_BLOCK = 0
LUMP_DIAGONAL = False
PROJECT_RHS = False


@register_solver(SOLVER_SIMPLE)
@register_solver(SOLVER_PISO)
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
            self.hydrostatic_pressure = setup_hydrostatic_pressure(
                simulation, needs_initial_value=True
            )
        ph_every_timestep = 'p_hydrostatic' in sim.data

        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])

        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)

        # Get matrices
        Matrices = EQUATION_SUBTYPES[self.equation_subtype]
        matrices = Matrices(
            simulation,
            use_stress_divergence_form=self.use_stress_divergence_form,
            use_grad_p_form=self.use_grad_p_form,
            use_grad_q_form=self.use_grad_q_form,
            use_lagrange_multiplicator=self.use_lagrange_multiplicator,
            include_hydrostatic_pressure=ph_every_timestep,
            incompressibility_flux_type=self.incompressibility_flux_type,
            num_elements_in_block=self.num_elements_in_block,
            lump_diagonal=self.lump_diagonal,
        )
        self.matrices = matrices

        # Slope limiter for the momentum equation velocity components
        self.slope_limiter = SlopeLimiterVelocity(
            sim, sim.data['u'], 'u', vel_w=sim.data['u_conv']
        )
        self.using_limiter = self.slope_limiter.active

        # Projection for the velocity
        self.velocity_postprocessor = None
        if self.velocity_postprocessing == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(
                sim,
                sim.data['u'],
                incompressibility_flux_type=self.incompressibility_flux_type,
            )

        # Storage for preassembled matrices
        self.A = None
        self.A_tilde = None
        self.A_tilde_inv = None
        self.B = None
        self.C = None

        # Temporary matrices to store matrix matrix products
        self.mat_AinvB = None  # SIMPLE & PISO
        self.mat_CAinvB = None  # SIMPLE & PISO
        self.mat_AinvA = None  # PISO
        self.mat_CAinvA = None  # PISO

        # Store number of iterations
        self.niters_u = None
        self.niters_p = None

    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation

        # PISO mode switch
        self.solver_type = sim.input.get_value('solver/type', required_type='string')

        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(
            self.simulation, 'solver/u', default_parameters=SOLVER_U_OPTIONS
        )
        self.pressure_solver = linear_solver_from_input(
            self.simulation, 'solver/p', default_parameters=SOLVER_P_OPTIONS
        )

        # Get the class to be used for the equation system assembly
        self.equation_subtype = sim.input.get_value(
            'solver/equation_subtype', EQUATION_SUBTYPE, 'string'
        )
        verify_key(
            'equation sub-type',
            self.equation_subtype,
            EQUATION_SUBTYPES,
            'SIMPLE solver',
        )

        # Lagrange multiplicator or remove null space via PETSc
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value(
            'solver/use_lagrange_multiplicator', USE_LAGRANGE_MULTIPLICATOR, 'bool'
        )
        if self.use_lagrange_multiplicator:
            self.remove_null_space = False

        # No need for special treatment if the pressure is set via Dirichlet conditions somewhere
        # No need for any tricks if the pressure is set via Dirichlet conditions somewhere
        if sim.data['dirichlet_bcs'].get('p', []) or sim.data['outlet_bcs']:
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

        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = Vu_family == 'Discontinuous Lagrange'

        # Velocity post_processing
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing = sim.input.get_value(
            'solver/velocity_postprocessing', default_postprocessing, 'string'
        )
        self.project_rhs = sim.input.get_value(
            'solver/project_rhs', PROJECT_RHS, 'bool'
        )
        verify_key(
            'velocity post processing',
            self.velocity_postprocessing,
            ('none', BDM),
            'SIMPLE solver',
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
    def momentum_prediction(self):
        """
        Solve the momentum prediction equation
        """
        sim = self.simulation

        def mom_assemble():
            """
            Assemble the linear systems
            """
            p_star = sim.data['p']
            if self.solver_type == SOLVER_SIMPLE:
                alpha = sim.input.get_value('solver/relaxation_u', ALPHA_U, 'float')
            else:
                alpha = 1.0
            assert 0 < alpha
            relax = (1 - alpha) / alpha

            # Get coupled guess velocities
            uvw_star = sim.data['uvw_star']

            # Assemble the LHS matrices only the first inner iteration
            if self.inner_iteration == 1:
                A, A_tilde, A_tilde_inv, B, C = self.matrices.assemble_matrices()
                self.A = dolfin.as_backend_type(A)
                self.A_tilde = dolfin.as_backend_type(A_tilde)
                self.A_tilde_inv = dolfin.as_backend_type(A_tilde_inv)
                self.B = dolfin.as_backend_type(B)
                self.C = dolfin.as_backend_type(C)

            A = self.A
            A_tilde = self.A_tilde
            B = self.B
            D = dolfin.as_backend_type(self.matrices.assemble_D())
            self.D = D

            if alpha != 1.0:
                lhs = dolfin.as_backend_type(A.copy())
                lhs.axpy(relax, A_tilde, False)
                lhs.apply('insert')
                rhs = D
                rhs.axpy(-1.0, B * p_star.vector())
                rhs.axpy(relax, A_tilde * uvw_star.vector())
                rhs.apply('insert')

            else:
                lhs = A
                rhs = D - B * p_star.vector()

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
                self.rhs_postprocessor = VelocityBDMProjection(
                    sim,
                    self.rhs_tmp,
                    incompressibility_flux_type=self.incompressibility_flux_type,
                )

            self.assigner_split.assign(list(self.rhs_tmp), rhs)
            self.rhs_postprocessor.run()
            self.assigner_merge.assign(rhs, list(self.rhs_tmp))

        def mom_solve(lhs, rhs):
            """
            Solve the linear systems
            """
            u_star = sim.data['uvw_star']
            u_temp = sim.data['uvw_temp']
            u_temp.assign(u_star)
            self.niters_u = self.velocity_solver.inner_solve(
                lhs,
                u_star.vector(),
                rhs,
                in_iter=self.inner_iteration,
                co_iter=self.co_inner_iter,
            )

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
    def pressure_correction(self, piso_rhs=False):
        """
        Solve the Navier-Stokes equations on SIMPLE form
        (Semi-Implicit Method for Pressure-Linked Equations)
        """
        sim = self.simulation
        p_hat = sim.data['p_hat']
        if self.solver_type == SOLVER_SIMPLE:
            alpha = sim.input.get_value('solver/relaxation_p', ALPHA_P, 'float')
        else:
            alpha = 1.0

        # Compute the LHS = C⋅Ãinv⋅B
        if self.inner_iteration == 1:
            C, Ainv, B = self.C, self.A_tilde_inv, self.B
            self.mat_AinvB = matmul(Ainv, B, self.mat_AinvB)
            self.mat_CAinvB = matmul(C, self.mat_AinvB, self.mat_CAinvB)
            self.LHS_pressure = dolfin.as_backend_type(self.mat_CAinvB.copy())
        LHS = self.LHS_pressure

        # Compute the RHS
        if not piso_rhs:
            # Standard SIMPLE pressure correction
            # Compute the divergence of u* and the rest of the right hand side
            uvw_star = sim.data['uvw_star']
            RHS = self.matrices.assemble_E_star(uvw_star)
            self.niters_p = 0
        else:
            # PISO pressure correction (the second pressure correction)
            # Compute the RHS = - C⋅Ãinv⋅(Ãinv - A)⋅û
            C, Ainv, A = self.C, self.A_tilde_inv, self.A
            if self.inner_iteration == 1:
                self.mat_AinvA = matmul(Ainv, A, self.mat_AinvA)
                self.mat_CAinvA = matmul(C, self.mat_AinvA, self.mat_CAinvA)
            RHS = self.mat_CAinvA * self.minus_uvw_hat - C * self.minus_uvw_hat

        # Inform PETSc about the null space
        if self.remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p_hat.vector())
                null_vec[:] = 1
                null_vec *= 1 / null_vec.norm("l2")

                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])

            # Make sure the null space is set on the matrix
            if self.inner_iteration == 1:
                LHS.set_nullspace(self.pressure_null_space)

            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(RHS)

        # Solve for the new pressure correction
        self.niters_p += self.pressure_solver.inner_solve(
            LHS,
            p_hat.vector(),
            RHS,
            in_iter=self.inner_iteration,
            co_iter=self.co_inner_iter,
        )

        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p_hat.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1) * dx2)
            pavg = dolfin.assemble(p_hat * dx2) / vol
            p_hat.vector()[:] -= pavg

        # Calculate p = p* + α p^
        sim.data['p'].vector().axpy(alpha, p_hat.vector())
        sim.data['p'].vector().apply('insert')

        return p_hat.vector().norm('l2')

    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """
        uvw = self.simulation.data['uvw_star']
        p_hat = self.simulation.data['p_hat']
        minus_uvw_hat = self.mat_AinvB * p_hat.vector()
        uvw.vector().axpy(-1.0, minus_uvw_hat)
        uvw.vector().apply('insert')

        if self.solver_type == SOLVER_PISO:
            self.minus_uvw_hat = minus_uvw_hat

    @timeit
    def velocity_update_piso(self):
        """
        Compute the final PISO velocity u*** using the twice corrected
        PISO pressure p***
        """
        sim = self.simulation

        p_sss = sim.data['p']
        u_ss = sim.data['uvw_star']

        # Compute the second velocity hat which is û² = Â⁻¹ (d - B p*** - A u**)
        u_hat2 = self.A_tilde_inv * (
            self.D - self.B * p_sss.vector() - self.A * u_ss.vector()
        )

        # Update the velocity from u** to u***
        u_ss.vector().axpy(1.0, u_hat2)
        u_ss.vector().apply('insert')

        self.uvw_hat2 = u_hat2
        return u_hat2.norm('l2')

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

    @timeit.named('run SIMPLE/PISO solver')
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

        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            num_inner_iter = sim.input.get_value(
                'solver/num_inner_iter', MAX_INNER_ITER, 'int'
            )
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

            # Update the coupled version of the velocity field
            self.assigner_merge.assign(sim.data['uvw_star'], list(sim.data['u']))

            # Run inner iterations
            self.inner_iteration = 1
            while self.inner_iteration <= num_inner_iter:
                self.co_inner_iter = num_inner_iter - self.inner_iteration

                if self.solver_type == SOLVER_SIMPLE:
                    err_u = self.momentum_prediction()
                    err_p = self.pressure_correction()
                    self.velocity_update()

                elif self.solver_type == SOLVER_PISO:
                    if self.inner_iteration == 1:
                        self.momentum_prediction()
                        self.pressure_correction()
                        self.velocity_update()
                    else:
                        self.niters_u = self.niters_p = 0
                        self.minus_uvw_hat = -1.0 * self.uvw_hat2
                    err_p = self.pressure_correction(piso_rhs=True)
                    err_u = self.velocity_update_piso()

                sim.log.info(
                    '  %s iteration %3d - err u* %10.3e - err p %10.3e'
                    ' - Num Krylov iters - u %3d - p %3d'
                    % (
                        self.solver_type,
                        self.inner_iteration,
                        err_u,
                        err_p,
                        self.niters_u,
                        self.niters_p,
                    )
                )
                self.inner_iteration += 1

                if err_u < allowable_error_inner:
                    break

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
