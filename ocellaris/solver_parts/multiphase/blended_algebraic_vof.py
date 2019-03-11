# Copyright (C) 2014-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from dolfin import Function, Constant
from ocellaris.solver_parts import SlopeLimiter
from . import register_multi_phase_model, MultiPhaseModel
from ..convection import get_convection_scheme, StaticScheme, VelocityDGT0Projector
from ocellaris.utils import linear_solver_from_input
from .vof import VOFMixin
from .advection_equation import AdvectionEquation


# Default values, can be changed in the input file
CONVECTION_SCHEME = 'Upwind'
CONTINUOUS_FIELDS = False
CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION = False
FORCE_STATIC = False
FORCE_BOUNDED = False
FORCE_SHARP = False
PLOT_FIELDS = False
NUM_SUBCYCLES = 1


# Default values, can be changed in the input file
SOLVER_OPTIONS = {
    'use_ksp': True,
    'petsc_ksp_type': 'gmres',
    'petsc_pc_type': 'asm',
    'petsc_ksp_initial_guess_nonzero': True,
    'petsc_ksp_view': 'DISABLED',
    'inner_iter_rtol': [1e-10] * 3,
    'inner_iter_atol': [1e-15] * 3,
    'inner_iter_max_it': [1000] * 3,
}


@register_multi_phase_model('BlendedAlgebraicVOF')
class BlendedAlgebraicVofModel(VOFMixin, MultiPhaseModel):
    description = 'A blended algebraic VOF scheme implementing HRIC/CICSAM type schemes'

    def __init__(self, simulation):
        """
        A blended algebraic VOF scheme works by using a specific
        convection scheme in the advection of the colour function
        that ensures a sharp interface.

        * The convection scheme should be the name of a convection
          scheme that is tailored for advection of the colour
          function, i.e "HRIC", "MHRIC", "RHRIC" etc,
        * The velocity field should be divergence free

        The colour function is unity when rho=rho0 and nu=nu0 and
        zero when rho=rho1 and nu=nu1
        """
        self.simulation = simulation
        simulation.log.info('Creating blended VOF multiphase model')

        # Define function space and solution function
        V = simulation.data['Vc']
        self.degree = V.ufl_element().degree()
        simulation.data['c'] = Function(V)
        simulation.data['cp'] = Function(V)
        simulation.data['cpp'] = Function(V)

        # The projected density and viscosity functions for the new time step can be made continuous
        self.continuous_fields = simulation.input.get_value(
            'multiphase_solver/continuous_fields', CONTINUOUS_FIELDS, 'bool'
        )
        if self.continuous_fields:
            simulation.log.info('    Using continuous rho and nu fields')
            mesh = simulation.data['mesh']
            V_cont = dolfin.FunctionSpace(mesh, 'CG', self.degree + 1)
            self.continuous_c = dolfin.Function(V_cont)
            self.continuous_c_old = dolfin.Function(V_cont)
            self.continuous_c_oldold = dolfin.Function(V_cont)

        self.force_bounded = simulation.input.get_value(
            'multiphase_solver/force_bounded', FORCE_BOUNDED, 'bool'
        )
        self.force_sharp = simulation.input.get_value(
            'multiphase_solver/force_sharp', FORCE_SHARP, 'bool'
        )

        # Calculate mu from rho and nu (i.e mu is quadratic in c) or directly from c (linear in c)
        self.calculate_mu_directly_from_colour_function = simulation.input.get_value(
            'multiphase_solver/calculate_mu_directly_from_colour_function',
            CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION,
            'bool',
        )

        # Get the physical properties
        self.set_physical_properties(read_input=True)

        # The convection blending function that counteracts numerical diffusion
        scheme = simulation.input.get_value(
            'convection/c/convection_scheme', CONVECTION_SCHEME, 'string'
        )
        simulation.log.info('    Using convection scheme %s for the colour function' % scheme)
        scheme_class = get_convection_scheme(scheme)
        self.convection_scheme = scheme_class(simulation, 'c')
        self.need_gradient = scheme_class.need_alpha_gradient

        # Create the equations when the simulation starts
        simulation.hooks.add_pre_simulation_hook(
            self.on_simulation_start, 'BlendedAlgebraicVofModel setup equations'
        )

        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(
            self.update, 'BlendedAlgebraicVofModel - update colour field'
        )
        simulation.hooks.register_custom_hook_point('MultiPhaseModelUpdated')

        # Linear solver
        # This causes the MPI unit tests to fail in "random" places for some reason
        # Quick fix: lazy loading of the solver
        LAZY_LOAD_SOLVER = True
        if LAZY_LOAD_SOLVER:
            self.solver = None
        else:
            self.solver = linear_solver_from_input(
                self.simulation, 'solver/c', default_parameters=SOLVER_OPTIONS
            )

        # Subcycle the VOF calculation multiple times per Navier-Stokes time step
        self.num_subcycles = scheme = simulation.input.get_value(
            'multiphase_solver/num_subcycles', NUM_SUBCYCLES, 'int'
        )
        if self.num_subcycles < 1:
            self.num_subcycles = 1

        # Time stepping based on the subcycled values
        if self.num_subcycles == 1:
            self.cp = simulation.data['cp']
            self.cpp = simulation.data['cpp']
        else:
            self.cp = dolfin.Function(V)
            self.cpp = dolfin.Function(V)

        # Plot density and viscosity fields for visualization
        self.plot_fields = simulation.input.get_value(
            'multiphase_solver/plot_fields', PLOT_FIELDS, 'bool'
        )
        if self.plot_fields:
            V_plot = V if not self.continuous_fields else V_cont
            self.rho_for_plot = Function(V_plot)
            self.nu_for_plot = Function(V_plot)
            self.rho_for_plot.rename('rho', 'Density')
            self.nu_for_plot.rename('nu', 'Kinematic viscosity')
            simulation.io.add_extra_output_function(self.rho_for_plot)
            simulation.io.add_extra_output_function(self.nu_for_plot)

        # Slope limiter in case we are using DG1, not DG0
        self.slope_limiter = SlopeLimiter(simulation, 'c', simulation.data['c'])
        simulation.log.info('    Using slope limiter: %s' % self.slope_limiter.limiter_method)
        self.is_first_timestep = True

    def on_simulation_start(self):
        """
        This runs when the simulation starts. It does not run in __init__
        since the solver needs the density and viscosity we define, and
        we need the velocity that is defined by the solver
        """
        sim = self.simulation
        beta = self.convection_scheme.blending_function

        # The time step (real value to be supplied later)
        self.dt = Constant(sim.dt / self.num_subcycles)

        # Setup the equation to solve
        c = sim.data['c']
        cp = self.cp
        cpp = self.cpp
        dirichlet_bcs = sim.data['dirichlet_bcs'].get('c', [])

        # Use backward Euler (BDF1) for timestep 1
        self.time_coeffs = Constant([1, -1, 0])

        if dolfin.norm(cpp.vector()) > 0 and self.num_subcycles == 1:
            # Use BDF2 from the start
            self.time_coeffs.assign(Constant([3 / 2, -2, 1 / 2]))
            sim.log.info('Using second order timestepping from the start in BlendedAlgebraicVOF')

        # Make sure the convection scheme has something useful in the first iteration
        c.assign(sim.data['cp'])

        if self.num_subcycles > 1:
            cp.assign(sim.data['cp'])

        # Plot density and viscosity
        self.update_plot_fields()

        # Define equation for advection of the colour function
        #    ∂c/∂t +  ∇⋅(c u) = 0
        Vc = sim.data['Vc']
        project_dgt0 = sim.input.get_value('multiphase_solver/project_uconv_dgt0', True, 'bool')
        if self.degree == 0 and project_dgt0:
            self.vel_dgt0_projector = VelocityDGT0Projector(sim, sim.data['u_conv'])
            self.u_conv = self.vel_dgt0_projector.velocity
        else:
            self.u_conv = sim.data['u_conv']
        forcing_zones = sim.data['forcing_zones'].get('c', [])
        self.eq = AdvectionEquation(
            sim,
            Vc,
            cp,
            cpp,
            self.u_conv,
            beta,
            time_coeffs=self.time_coeffs,
            dirichlet_bcs=dirichlet_bcs,
            forcing_zones=forcing_zones,
            dt=self.dt,
        )

        if self.need_gradient:
            # Reconstruct the gradient from the colour function DG0 field
            self.convection_scheme.initialize_gradient()

        # Notify listeners that the initial values are available
        sim.hooks.run_custom_hook('MultiPhaseModelUpdated')

    def get_colour_function(self, k):
        """
        Return the colour function on timestep t^{n+k}
        """
        if k == 0:
            if self.continuous_fields:
                c = self.continuous_c
            else:
                c = self.simulation.data['c']
        elif k == -1:
            if self.continuous_fields:
                c = self.continuous_c_old
            else:
                c = self.simulation.data['cp']
        elif k == -2:
            if self.continuous_fields:
                c = self.continuous_c_oldold
            else:
                c = self.simulation.data['cpp']

        if self.force_bounded:
            c = dolfin.max_value(dolfin.min_value(c, Constant(1.0)), Constant(0.0))

        if self.force_sharp:
            c = dolfin.conditional(dolfin.ge(c, 0.5), Constant(1.0), Constant(0.0))

        return c

    def update_plot_fields(self):
        """
        These fields are only needed to visualise the rho and nu fields
        in xdmf format for Paraview or similar
        """
        if not self.plot_fields:
            return
        V = self.rho_for_plot.function_space()
        dolfin.project(self.get_density(0), V, function=self.rho_for_plot)
        dolfin.project(self.get_laminar_kinematic_viscosity(0), V, function=self.nu_for_plot)

    def update(self, timestep_number, t, dt):
        """
        Update the VOF field by advecting it for a time dt
        using the given divergence free velocity field
        """
        timer = dolfin.Timer('Ocellaris update VOF')
        sim = self.simulation

        # Get the functions
        c = sim.data['c']
        cp = sim.data['cp']
        cpp = sim.data['cpp']

        # Stop early if the free surface is forced to stay still
        force_static = sim.input.get_value('multiphase_solver/force_static', FORCE_STATIC, 'bool')
        if force_static:
            c.assign(cp)
            cpp.assign(cp)
            timer.stop()  # Stop timer before hook
            sim.hooks.run_custom_hook('MultiPhaseModelUpdated')
            self.is_first_timestep = False
            return

        if timestep_number != 1:
            # Update the previous values
            cpp.assign(cp)
            cp.assign(c)

            if self.degree == 0:
                self.vel_dgt0_projector.update()

        # Reconstruct the gradients
        if self.need_gradient:
            self.convection_scheme.gradient_reconstructor.reconstruct()

        # Update the convection blending factors
        is_static = isinstance(self.convection_scheme, StaticScheme)
        if not is_static:
            self.convection_scheme.update(dt / self.num_subcycles, self.u_conv)

        # Update global bounds in slope limiter
        if self.is_first_timestep:
            lo, hi = self.slope_limiter.set_global_bounds(lo=0.0, hi=1.0)
            if self.slope_limiter.has_global_bounds:
                sim.log.info(
                    'Setting global bounds [%r, %r] in BlendedAlgebraicVofModel' % (lo, hi)
                )

        # Solve the advection equations for the colour field
        if timestep_number == 1 or is_static:
            c.assign(cp)
        else:
            if self.solver is None:
                sim.log.info('Creating colour function solver', flush=True)
                self.solver = linear_solver_from_input(
                    self.simulation, 'solver/c', default_parameters=SOLVER_OPTIONS
                )

            # Solve the advection equation
            A = self.eq.assemble_lhs()
            for _ in range(self.num_subcycles):
                b = self.eq.assemble_rhs()
                self.solver.inner_solve(A, c.vector(), b, 1, 0)
                self.slope_limiter.run()
                if self.num_subcycles > 1:
                    self.cpp.assign(self.cp)
                    self.cp.assign(c)

        # Optionally use a continuous predicted colour field
        if self.continuous_fields:
            Vcg = self.continuous_c.function_space()
            dolfin.project(c, Vcg, function=self.continuous_c)
            dolfin.project(cp, Vcg, function=self.continuous_c_old)
            dolfin.project(cpp, Vcg, function=self.continuous_c_oldold)

        # Report properties of the colour field
        sim.reporting.report_timestep_value('min(c)', c.vector().min())
        sim.reporting.report_timestep_value('max(c)', c.vector().max())

        # The next update should use the dt from this time step of the
        # main Navier-Stoke solver. The update just computed above uses
        # data from the previous Navier-Stokes solve with the previous dt
        self.dt.assign(dt / self.num_subcycles)

        if dt != sim.dt_prev:
            # Temporary switch to first order timestepping for the next
            # time step. This code is run before the Navier-Stokes solver
            # in each time step
            sim.log.info('VOF solver is first order this time step due to change in dt')
            self.time_coeffs.assign(Constant([1.0, -1.0, 0.0]))
        else:
            # Use second order backward time difference next time step
            self.time_coeffs.assign(Constant([3 / 2, -2.0, 1 / 2]))

        self.update_plot_fields()
        timer.stop()  # Stop timer before hook
        sim.hooks.run_custom_hook('MultiPhaseModelUpdated')
        self.is_first_timestep = False
