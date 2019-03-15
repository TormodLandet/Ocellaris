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
)
from . import Solver, register_solver, BDM
from ..solver_parts import (
    VelocityBDMProjection,
    setup_hydrostatic_pressure,
    SlopeLimiterVelocity,
    before_simulation,
    after_timestep,
)
from .ipcs_equations import EQUATION_SUBTYPES


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
HYDROSTATIC_PRESSURE_CALCULATION_EVERY_TIMESTEP = False
INCOMPRESSIBILITY_FLUX_TYPE = 'central'


@register_solver('IPCS-D')
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
        (
            MomentumPredictionEquation,
            PressureCorrectionEquation,
            VelocityUpdateEquation,
        ) = EQUATION_SUBTYPES[self.equation_subtype]

        # Define the momentum prediction equations
        self.eqs_mom_pred = MomentumPredictionEquation(
            simulation,
            use_stress_divergence_form=self.use_stress_divergence_form,
            use_grad_p_form=self.use_grad_p_form,
            include_hydrostatic_pressure=self.hydrostatic_pressure.every_timestep,
        )

        # Define the pressure correction equation
        self.eq_pressure = PressureCorrectionEquation(
            simulation,
            use_lagrange_multiplicator=self.use_lagrange_multiplicator,
            incompressibility_flux_type=self.incompressibility_flux_type,
        )

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
            self.velocity_postprocessor = VelocityBDMProjection(
                sim, sim.data['u'], incompressibility_flux_type=self.incompressibility_flux_type
            )

        # Storage for preassembled matrices
        self.Au = None
        self.Ap = None
        self.Au_upd = None

        # Store number of iterations
        self.niters_u = None
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
        self.vel_is_discontinuous = Vu_family == 'Discontinuous Lagrange'

        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(
            self.simulation, 'solver/u', default_parameters=SOLVER_U_OPTIONS
        )
        self.pressure_solver = linear_solver_from_input(
            self.simulation, 'solver/p', default_parameters=SOLVER_P_OPTIONS
        )

        # Velocity update can be performed with local solver for DG velocities
        self.use_local_solver_for_update = sim.input.get_value(
            'solver/u_upd_local', self.vel_is_discontinuous, 'bool'
        )
        if self.use_local_solver_for_update:
            self.u_upd_solver = None  # Will be set when LHS is ready
        else:
            self.u_upd_solver = linear_solver_from_input(
                self.simulation, 'solver/u_upd', default_parameters=SOLVER_U_OPTIONS
            )

        # Get the class to be used for the equation system assembly
        self.equation_subtype = sim.input.get_value(
            'solver/equation_subtype', EQUATION_SUBTYPE, 'string'
        )
        verify_key('equation sub-type', self.equation_subtype, EQUATION_SUBTYPES, 'ipcs solver')

        # Lagrange multiplicator or remove null space via PETSc
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value(
            'solver/use_lagrange_multiplicator', USE_LAGRANGE_MULTIPLICATOR, 'bool'
        )
        has_dirichlet = self.simulation.data['dirichlet_bcs'].get('p', []) or sim.data['outlet_bcs']
        if self.use_lagrange_multiplicator or has_dirichlet:
            self.remove_null_space = False

        # No need for special treatment if the pressure is set via Dirichlet conditions somewhere
        if has_dirichlet:
            self.use_lagrange_multiplicator = False
            self.remove_null_space = False

        # Control the form of the governing equations
        self.use_stress_divergence_form = sim.input.get_value(
            'solver/use_stress_divergence_form', USE_STRESS_DIVERGENCE, 'bool'
        )
        self.use_grad_p_form = sim.input.get_value(
            'solver/use_grad_p_form', USE_GRAD_P_FORM, 'bool'
        )
        self.incompressibility_flux_type = sim.input.get_value(
            'solver/incompressibility_flux_type', INCOMPRESSIBILITY_FLUX_TYPE, 'string'
        )

        # Velocity post_processing
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing = sim.input.get_value(
            'solver/velocity_postprocessing', default_postprocessing, 'string'
        )
        verify_key(
            'velocity post processing', self.velocity_postprocessing, (None, BDM), 'ipcs solver'
        )

        # Quasi-steady simulation input
        self.steady_velocity_eps = sim.input.get_value(
            'solver/steady_velocity_stopping_criterion', None, 'float'
        )
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
    def momentum_prediction(self):
        """
        Solve the momentum prediction equation
        """
        sim = self.simulation
        eq = self.eqs_mom_pred

        # Collect previous velocity components in coupled function
        uvw_star = sim.data['uvw_star']
        uvw_temp = sim.data['uvw_temp']
        uvw_temp.assign(uvw_star)

        if self.inner_iteration == 1:
            # Assemble the A matrix only the first inner iteration
            self.Au = dolfin.as_backend_type(eq.assemble_lhs())

        A = self.Au
        b = dolfin.as_backend_type(eq.assemble_rhs())
        self.assigner_merge.assign(uvw_star, list(sim.data['u']))
        self.niters_u = self.velocity_solver.inner_solve(
            A, uvw_star.vector(), b, in_iter=self.inner_iteration, co_iter=self.co_inner_iter
        )
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
        if self.inner_iteration == 1:
            self.Ap = dolfin.as_backend_type(self.eq_pressure.assemble_lhs())
        A = self.Ap
        b = dolfin.as_backend_type(self.eq_pressure.assemble_rhs())

        # Inform PETSc about the null space
        if self.remove_null_space:
            if self.pressure_null_space is None:
                # Create vector that spans the null space
                null_vec = dolfin.Vector(p.vector())
                null_vec[:] = 1
                null_vec *= 1 / null_vec.norm("l2")

                # Create null space basis object
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])

            # Make sure the null space is set on the matrix
            if self.inner_iteration == 1:
                A.set_nullspace(self.pressure_null_space)

            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(b)

        # Solve for the new pressure correction
        self.niters_p = self.pressure_solver.inner_solve(
            A, p.vector(), b, in_iter=self.inner_iteration, co_iter=self.co_inner_iter
        )

        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1) * dx2)
            pavg = dolfin.assemble(p * dx2) / vol
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
                eq = self.eqs_vel_upd[d]

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

            # Run inner iterations
            self.inner_iteration = 1
            while self.inner_iteration <= num_inner_iter:
                self.co_inner_iter = num_inner_iter - self.inner_iteration
                err_u = self.momentum_prediction()
                err_p = self.pressure_correction()
                self.velocity_update()
                sim.log.info(
                    '  IPCS-D iteration %3d - err u* %10.3e - err p %10.3e'
                    ' - Num Krylov iters - u %3d - p %3d'
                    % (self.inner_iteration, err_u, err_p, self.niters_u, self.niters_p)
                )
                self.inner_iteration += 1

                if err_u < allowable_error_inner:
                    break

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
