# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from dolfin import dot, grad, dx, sqrt, Form, assemble, Timer
from dolfin import TrialFunction, TestFunction, FunctionSpace, Function
from dolfin import VectorSpaceBasis, as_backend_type
from ocellaris.utils import linear_solver_from_input


DEFAULT_SOLVER_CONFIGURATION = {
    'use_ksp': True,
    # Set PETSc solve type (conjugate gradient) and preconditioner
    # (algebraic multigrid)
    'petsc_ksp_type': 'cg',
    'petsc_pc_type': 'gamg',
    # Since we have a singular problem, use SVD solver on the multigrid
    # 'coarse grid'
    'petsc_mg_coarse_ksp_type': 'preonly',
    'petsc_mg_coarse_pc_type': 'svd',
    # Set the solver tolerance and allow starting from previous solution
    'petsc_ksp_rtol': 1e-10,
    'petsc_ksp_initial_guess_nonzero': True,
    # Verbosity
    'petsc_ksp_view': 'DISABLED',
    'petsc_ksp_monitor': 'DISABLED',
}


class HydrostaticPressure:
    def __init__(self, simulation, every_timestep):
        """
        Calculate the hydrostatic pressure
        """
        self.simulation = simulation
        self.active = True
        self.every_timestep = every_timestep
        rho = simulation.data['rho']
        g = simulation.data['g']
        ph = simulation.data['p_hydrostatic']

        # Define the weak form
        Vp = ph.function_space()
        p = TrialFunction(Vp)
        q = TestFunction(Vp)
        a = dot(grad(p), grad(q)) * dx
        L = rho * dot(g, grad(q)) * dx

        self.func = ph
        self.tensor_lhs = assemble(a)
        self.form_rhs = Form(L)
        self.null_space = None
        self.solver = linear_solver_from_input(
            simulation,
            'solver/p_hydrostatic',
            default_parameters=DEFAULT_SOLVER_CONFIGURATION,
        )

    def update(self):
        if not self.active:
            return

        with Timer('Ocellaris update hydrostatic pressure'):
            A = self.tensor_lhs
            b = assemble(self.form_rhs)

            if self.null_space is None:
                # Create vector that spans the null space, and normalize
                null_space_vector = b.copy()
                null_space_vector[:] = sqrt(1.0 / null_space_vector.size())

                # Create null space basis object and attach to PETSc matrix
                self.null_space = VectorSpaceBasis([null_space_vector])
                as_backend_type(A).set_nullspace(self.null_space)

            self.null_space.orthogonalize(b)
            self.solver.solve(A, self.func.vector(), b)

        if not self.every_timestep:
            # Give initial values for p, but do not continuously compute p_hydrostatic
            sim = self.simulation
            p = sim.data['p']
            if p.vector().max() == p.vector().min() == 0.0:
                sim.log.info(
                    'Initial pressure field is identically zero, initializing to hydrostatic'
                )
                p.interpolate(self.func)

            # Disable further hydrostatic pressure calculations
            self.func.vector().zero()
            del sim.data['p_hydrostatic']
            del self.func
            self.active = False


class NoHydrostaticPressure:
    every_timestep = False

    def update(self):
        pass


def setup_hydrostatic_pressure(
    simulation, needs_initial_value, default_every_timestep=False
):
    """
    We can calculate the hydrostatic pressure as its own pressure field every
    time step such that the we only solves for the dynamic pressure. For
    segregated solvers we
    """
    # No hydrostatic pressure field if there is no gravity
    has_gravity = any(gi != 0 for gi in simulation.data['g'].values())
    if not has_gravity:
        return NoHydrostaticPressure()

    # We only calculate the hydrostatic pressure every time step if asked
    ph_every_timestep = simulation.input.get_value(
        'solver/hydrostatic_pressure_calculation_every_timestep',
        default_every_timestep,
        required_type='bool',
    )

    if not (needs_initial_value or ph_every_timestep):
        return NoHydrostaticPressure()

    # Hydrostatic pressure is always CG
    simulation.log.info('    Setting up hydrostatic pressure calculations')
    Vp = simulation.data['Vp']
    Pp = Vp.ufl_element().degree()
    Vph = FunctionSpace(simulation.data['mesh'], 'CG', Pp)
    simulation.data['p_hydrostatic'] = Function(Vph)

    # Helper class to calculate the hydrostatic pressure distribution
    hydrostatic_pressure = HydrostaticPressure(simulation, ph_every_timestep)

    return hydrostatic_pressure
