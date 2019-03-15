# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from dolfin import Constant
from ocellaris.utils import (
    ocellaris_error,
    timeit,
    linear_solver_from_input,
    create_vector_functions,
    shift_fields,
    velocity_change,
)
from ocellaris.solver_parts import (
    setup_hydrostatic_pressure,
    VelocityBDMProjection,
    SlopeLimiterVelocity,
    before_simulation,
    after_timestep,
)
from . import Solver, register_solver, BDM, UPWIND
from .coupled_equations import EQUATION_SUBTYPES


# Default values, can be changed in the input file
LU_SOLVER_1CPU = 'umfpack'
LU_SOLVER_NCPU = 'mumps'
LU_PARAMETERS = {}

# Default values, can be changed in the input file
EQUATION_SUBTYPE = 'Default'
USE_STRESS_DIVERGENCE = False
USE_LAGRANGE_MULTIPLICATOR = False
FIX_PRESSURE_DOF = True
USE_GRAD_P_FORM = False
USE_GRAD_Q_FORM = True
PRESSURE_CONTINUITY_FACTOR = 0
VELOCITY_CONTINUITY_FACTOR_D12 = 0
INCOMPRESSIBILITY_FLUX_TYPE = 'central'


@register_solver('Coupled')
class SolverCoupled(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes coupled solver based on DG-SIP

        :type simulation: ocellaris.Simulation
        """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()
        self.hydrostatic_pressure = setup_hydrostatic_pressure(
            simulation, needs_initial_value=False
        )

        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])

        # Solver control parameters
        sim.data['dt'] = Constant(simulation.dt)

        # Equation class
        CoupledEquations = EQUATION_SUBTYPES[self.equation_subtype]

        # Get the BCs for the coupled function space
        self.dirichlet_bcs = self.coupled_boundary_conditions(CoupledEquations.use_strong_bcs)

        # Create equation
        self.eqs = CoupledEquations(
            simulation,
            flux_type=self.flux_type,
            use_stress_divergence_form=self.use_stress_divergence_form,
            use_grad_p_form=self.use_grad_p_form,
            use_grad_q_form=self.use_grad_q_form,
            use_lagrange_multiplicator=self.use_lagrange_multiplicator,
            pressure_continuity_factor=self.pressure_continuity_factor,
            velocity_continuity_factor_D12=self.velocity_continuity_factor_D12,
            include_hydrostatic_pressure=self.hydrostatic_pressure.every_timestep,
            incompressibility_flux_type=self.incompressibility_flux_type,
        )

        sim.log.info('    Using velocity postprocessor: %r' % self.velocity_postprocessing_method)
        self.velocity_postprocessor = None
        if self.velocity_postprocessing_method == BDM:
            D12 = self.velocity_continuity_factor_D12
            self.velocity_postprocessor = VelocityBDMProjection(
                sim,
                sim.data['u'],
                incompressibility_flux_type=self.incompressibility_flux_type,
                D12=D12,
            )

        # Velocity slope limiter
        self.using_limiter = False
        if self.vel_is_discontinuous:
            self.slope_limiter = SlopeLimiterVelocity(
                sim, sim.data['u'], 'u', vel_w=sim.data['u_conv']
            )
            self.using_limiter = self.slope_limiter.active

        if self.fix_pressure_dof:
            pdof = get_global_row_number(self.subspaces[-1])
            self.pressure_row_to_fix = numpy.array([pdof], dtype=numpy.intc)

        # Store number of iterations
        self.niters = None

    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation

        # Solver for the coupled system
        default_lu_solver = LU_SOLVER_1CPU if sim.ncpu == 1 else LU_SOLVER_NCPU
        self.coupled_solver = linear_solver_from_input(
            sim, 'solver/coupled', 'lu', None, default_lu_solver, LU_PARAMETERS
        )

        # Get the class to be used for the equation system assembly
        self.equation_subtype = sim.input.get_value(
            'solver/equation_subtype', EQUATION_SUBTYPE, 'string'
        )
        if self.equation_subtype not in EQUATION_SUBTYPES:
            available_methods = '\n'.join(' - %s' % m for m in EQUATION_SUBTYPES)
            ocellaris_error(
                'Unknown equation sub-type',
                'Equation sub-type %s not available for coupled solver, please use one of:\n%s'
                % (self.equation_subtype, available_methods),
            )

        # Give warning if using iterative solver
        if isinstance(self.coupled_solver, dolfin.PETScKrylovSolver):
            sim.log.warning(
                'WARNING: Using a Krylov solver for the coupled NS equations is not a good idea'
            )
        else:
            # Removed in DOLFIN 2018.1, causes segfault now (Nov 2018)
            # self.coupled_solver.set_parameter('same_nonzero_pattern', True)
            pass

        # Lagrange multiplicator or remove null space via PETSc or just normalize after solving
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value(
            'solver/use_lagrange_multiplicator', USE_LAGRANGE_MULTIPLICATOR, 'bool'
        )
        self.fix_pressure_dof = sim.input.get_value(
            'solver/fix_pressure_dof', FIX_PRESSURE_DOF, 'bool'
        )
        if self.use_lagrange_multiplicator or self.fix_pressure_dof:
            self.remove_null_space = False

        # Check if the solver supports removing null spaces, otherwise we can enable a hack
        # that works better than nothing, but is most definitely not preferred
        self.normalize_pressure = False
        does_not_support_null_space = ('mumps',)
        if (
            self.remove_null_space
            and self.coupled_solver.created_with_lu_method in does_not_support_null_space
        ):
            self.normalize_pressure = True
            self.remove_null_space = False

        # No need for any tricks if the pressure is set via Dirichlet conditions somewhere
        if sim.data['dirichlet_bcs'].get('p', []) or sim.data['outlet_bcs']:
            self.remove_null_space = False
            self.use_lagrange_multiplicator = False
            self.fix_pressure_dof = False

        # Control the form of the governing equations
        self.flux_type = sim.input.get_value('convection/u/flux_type', UPWIND, 'string')
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
        self.pressure_continuity_factor = sim.input.get_value(
            'solver/pressure_continuity_factor', PRESSURE_CONTINUITY_FACTOR, 'float'
        )
        self.velocity_continuity_factor_D12 = sim.input.get_value(
            'solver/velocity_continuity_factor_D12', VELOCITY_CONTINUITY_FACTOR_D12, 'float'
        )

        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = Vu_family == 'Discontinuous Lagrange'

        # Local DG velocity postprocessing
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing_method = sim.input.get_value(
            'solver/velocity_postprocessing', default_postprocessing, 'string'
        )

    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation

        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        cd = sim.data['constrained_domain']

        # Create coupled mixed function space and mixed function to hold results
        func_spaces = [Vu] * sim.ndim + [Vp]
        self.subspace_names = ['u%d' % d for d in range(sim.ndim)] + ['p']

        if self.use_lagrange_multiplicator:
            Vl = dolfin.FunctionSpace(sim.data['mesh'], "R", 0, constrained_domain=cd)
            sim.data['l'] = dolfin.Function(Vl)
            func_spaces.append(Vl)
            self.subspace_names.append('l')

        e_mixed = dolfin.MixedElement([fs.ufl_element() for fs in func_spaces])
        Vcoupled = dolfin.FunctionSpace(sim.data['mesh'], e_mixed)
        sim.data['Vcoupled'] = Vcoupled
        sim.ndofs += Vcoupled.dim()

        Nspace = len(func_spaces)
        self.subspaces = [Vcoupled.sub(i) for i in range(Nspace)]
        sim.data['coupled'] = self.coupled_func = dolfin.Function(Vcoupled)
        self.assigner = dolfin.FunctionAssigner(func_spaces, Vcoupled)

        # Create segregated functions on component and vector form
        create_vector_functions(sim, 'u', 'u%d', Vu)
        create_vector_functions(sim, 'up', 'up%d', Vu)
        create_vector_functions(sim, 'upp', 'upp%d', Vu)
        create_vector_functions(sim, 'u_conv', 'u_conv%d', Vu)
        create_vector_functions(sim, 'up_conv', 'up_conv%d', Vu)
        create_vector_functions(sim, 'upp_conv', 'upp_conv%d', Vu)
        create_vector_functions(sim, 'u_unlim', 'u_unlim%d', Vu)
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['ui_tmp'] = dolfin.Function(Vu)

    def coupled_boundary_conditions(self, use_strong_bcs):
        """
        Convert boundary conditions from segregated to coupled function spaces
        """
        coupled_dirichlet_bcs = []
        for i, name in enumerate(self.subspace_names):
            if not use_strong_bcs and name.startswith('u'):
                # Use weak BCs if the velocity is DG
                continue

            V = self.subspaces[i]
            bcs = self.simulation.data['dirichlet_bcs'].get(name, [])
            for bc in bcs:
                bc_new = bc.copy_and_change_function_space(V)
                coupled_dirichlet_bcs.append(bc_new)

        return coupled_dirichlet_bcs

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

    @timeit
    def solve_coupled(self):
        """
        Solve the coupled equations
        """
        # Assemble the equation system
        A = self.eqs.assemble_lhs()
        b = self.eqs.assemble_rhs()

        # Apply strong boundary conditions (this list is empty for DG)
        for dbc in self.dirichlet_bcs:
            dbc.apply(A, b)

        if self.fix_pressure_dof:
            A.ident(self.pressure_row_to_fix)
        elif self.remove_null_space:
            if self.pressure_null_space is None:
                # Create null space vector in Vp Space
                null_func = dolfin.Function(self.simulation.data['Vp'])
                null_vec = null_func.vector()
                null_vec[:] = 1
                null_vec *= 1 / null_vec.norm("l2")

                # Convert null space vector to coupled space
                null_func2 = dolfin.Function(self.simulation.data['Vcoupled'])
                ndim = self.simulation.ndim
                fa = dolfin.FunctionAssigner(self.subspaces[ndim], self.simulation.data['Vp'])
                fa.assign(null_func2.sub(ndim), null_func)

                # Create the null space basis
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_func2.vector()])

            # Make sure the null space is set on the matrix
            dolfin.as_backend_type(A).set_nullspace(self.pressure_null_space)

            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(b)

        # Solve the equation system
        self.simulation.hooks.matrix_ready('Coupled', A, b)
        self.coupled_solver.solve(A, self.coupled_func.vector(), b)

        # Assign into the regular (split) functions from the coupled function
        funcs = [self.simulation.data[name] for name in self.subspace_names]
        self.assigner.assign(funcs, self.coupled_func)
        for func in funcs:
            func.vector().apply('insert')  # dolfin bug #587

        # Some solvers cannot remove the null space, so we just normalize the pressure instead.
        # If we remove the null space of the matrix system this will not be the exact same as
        # removing the proper null space of the equation, so we also fix this here
        if self.normalize_pressure or self.remove_null_space or self.fix_pressure_dof:
            p = self.simulation.data['p']
            dx2 = dolfin.dx(domain=p.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1) * dx2)
            # Perform correction multiple times due to round-of error. The first correction
            # can be i.e 1e14 while the next correction is around unity
            pavg = 1e10
            while abs(pavg) > 1000:
                pavg = dolfin.assemble(p * dx2) / vol
                p.vector()[:] -= pavg

    @timeit.named('run coupled solver')
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation
        sim.hooks.simulation_started()

        # Setup timestepping and initial convecting velocity
        force_steady = sim.input.get_value('solver/force_steady', False, 'bool')
        before_simulation(sim, force_steady)

        # Time loop
        t = sim.time
        it = sim.timestep
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            steady_eps = sim.input.get_value(
                'solver/steady_velocity_stopping_criterion', -1, 'float'
            )
            force_steady = sim.input.get_value('solver/force_steady', False, 'bool')

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

            # Solve for the new time step
            self.solve_coupled()

            # Postprocess the solution velocity field
            self.postprocess_velocity()

            # Set u_conv equal to u
            shift_fields(sim, ['u%d', 'u_conv%d'])

            # Slope limit the velocities
            if self.using_limiter:
                change_lim = self.slope_limit_velocities()
                sim.reporting.report_timestep_value('ulim_diff', change_lim)

            # Move u -> up, up -> upp and prepare for the next time step
            vel_diff = after_timestep(sim, steady_eps >= 0, force_steady)

            # Stop steady state simulation if convergence has been reached
            if steady_eps >= 0:
                vel_diff = dolfin.MPI.max(dolfin.MPI.comm_world, float(vel_diff))
                sim.reporting.report_timestep_value('max(ui_new-ui_prev)', vel_diff)
                if vel_diff < steady_eps:
                    sim.log.info('Stopping simulation, steady state achieved')
                    sim.input.set_value('time/tmax', t)

            # Postprocess this time step
            sim.hooks.end_timestep()

        # We are done
        sim.hooks.simulation_ended(success=True)


def get_global_row_number(V):
    """
    Get the lowest global matrix row number belonging to the
    function space V and local cell 0
    If V is a not a subspace then 0 will normally be returned
    (since dof 0 will typically belong to cell 0)
    """
    dm = V.dofmap()
    dof = dm.cell_dofs(0).min()
    gdof = dm.local_to_global_index(dof)
    return dolfin.MPI.min(dolfin.MPI.comm_world, int(gdof))
