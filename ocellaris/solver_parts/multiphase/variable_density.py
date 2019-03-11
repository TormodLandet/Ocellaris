# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from dolfin import Function, Constant
from . import register_multi_phase_model, MultiPhaseModel
from ocellaris.solver_parts import SlopeLimiter, RungeKuttaDGTimestepping
from ocellaris.utils import (
    linear_solver_from_input,
    ocellaris_interpolate,
    ocellaris_error,
)
from .advection_equation import AdvectionEquation


# Default values, can be changed in the input file
SOLVER = 'gmres'
PRECONDITIONER = 'default'
KRYLOV_PARAMETERS = {
    'nonzero_initial_guess': True,
    'relative_tolerance': 1e-15,
    'absolute_tolerance': 1e-15,
}
FORCE_STATIC = False


@register_multi_phase_model('VariableDensity')
class VariableDensityModel(MultiPhaseModel):
    description = 'A variable density DG transport equation model'

    def __init__(self, simulation):
        """
        A variable density scheme solving D/Dt(rho)=0 with a
        constant kinematic vosicsity nu and a variable dynamic
        visocisty mu (depending on rho)
        """
        self.simulation = simulation
        simulation.log.info('Creating variable density multiphase model')

        # Define function space and solution function
        V = simulation.data['Vrho']
        self.rho = simulation.data['rho'] = Function(V)
        self.rho_p = simulation.data['rho_p'] = Function(V)
        self.rho_pp = simulation.data['rho_pp'] = Function(V)

        # Get the physical properties
        inp = simulation.input
        self.rho_min = inp.get_value(
            'physical_properties/rho_min', required_type='float'
        )
        self.rho_max = inp.get_value(
            'physical_properties/rho_max', required_type='float'
        )
        self.nu = inp.get_value('physical_properties/nu', required_type='float')

        # Create the equations when the simulation starts
        simulation.hooks.add_pre_simulation_hook(
            self.on_simulation_start, 'VariableDensityModel setup equations'
        )

        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(
            self.update, 'VariableDensityModel - update density field'
        )
        simulation.hooks.register_custom_hook_point('MultiPhaseModelUpdated')

        self.use_analytical_solution = inp.get_value(
            'multiphase_solver/analytical_solution', False, 'bool'
        )
        self.use_rk_method = inp.get_value(
            'multiphase_solver/explicit_rk_method', False, 'bool'
        )
        self.is_first_timestep = True

    @classmethod
    def create_function_space(cls, simulation):
        mesh = simulation.data['mesh']
        cd = simulation.data['constrained_domain']
        Vr_name = simulation.input.get_value(
            'multiphase_solver/function_space_rho', 'Discontinuous Lagrange', 'string'
        )
        Pr = simulation.input.get_value(
            'multiphase_solver/polynomial_degree_rho', 1, 'int'
        )
        Vrho = dolfin.FunctionSpace(mesh, Vr_name, Pr, constrained_domain=cd)
        simulation.data['Vrho'] = Vrho
        simulation.ndofs += Vrho.dim()

    def on_simulation_start(self):
        """
        This runs when the simulation starts. It does not run in __init__
        since the N-S solver needs the density and viscosity we define, and
        we need the velocity that is defined by the solver
        """
        if self.use_analytical_solution:
            return

        sim = self.simulation
        dirichlet_bcs = sim.data['dirichlet_bcs'].get('rho', [])

        if self.use_rk_method:
            V = self.simulation.data['Vrho']
            if not V.ufl_element().family() == 'Discontinuous Lagrange':
                ocellaris_error(
                    'VariableDensity timestepping error',
                    'Can only use explicit SSP Runge-Kutta method with DG space for rho',
                )

            Vu = sim.data['Vu']
            u_conv, self.funcs_to_extrapolate = [], []
            for d in range(sim.ndim):
                ux = Function(Vu)
                up = sim.data['up_conv%d' % d]
                upp = sim.data['upp_conv%d' % d]
                self.funcs_to_extrapolate.append((ux, up, upp))
                u_conv.append(ux)
            u_conv = dolfin.as_vector(u_conv)

            from dolfin import dot, div, jump, dS

            mesh = self.simulation.data['mesh']

            re = self.rho_explicit = dolfin.Function(V)
            c, d = dolfin.TrialFunction(V), dolfin.TestFunction(V)
            n = dolfin.FacetNormal(mesh)
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n))) / 2
            dx = dolfin.dx(domain=mesh)

            eq = c * d * dx

            # Convection integrated by parts two times to bring back the original
            # div form (this means we must subtract and add all fluxes)
            eq += div(re * u_conv) * d * dx

            # Replace downwind flux with upwind flux on downwind internal facets
            eq -= jump(w_nD * d) * jump(re) * dS

            # Replace downwind flux with upwind BC flux on downwind external facets
            for dbc in dirichlet_bcs:
                # Subtract the "normal" downwind flux
                eq -= w_nD * re * d * dbc.ds()
                # Add the boundary value upwind flux
                eq += w_nD * dbc.func() * d * dbc.ds()

            a, L = dolfin.system(eq)
            self.rk = RungeKuttaDGTimestepping(
                self.simulation,
                a,
                L,
                self.rho,
                self.rho_explicit,
                'rho',
                order=None,
                explicit_funcs=self.funcs_to_extrapolate,
                bcs=dirichlet_bcs,
            )

        else:
            # Use backward Euler (BDF1) for timestep 1
            self.time_coeffs = Constant([1, -1, 0])

            if dolfin.norm(self.rho_pp.vector()) > 0:
                # Use BDF2 from the start
                self.time_coeffs.assign(Constant([3 / 2, -2, 1 / 2]))
                self.simulation.log.info(
                    'Using second order timestepping from the start in VariableDensity'
                )

            # Define equation for advection of the density
            #    ∂ρ/∂t +  ∇⋅(ρ u) = 0
            beta = None
            u_conv = sim.data['u_conv']
            forcing_zones = sim.data['forcing_zones'].get('rho', [])
            self.eq = AdvectionEquation(
                sim,
                sim.data['Vrho'],
                self.rho_p,
                self.rho_pp,
                u_conv,
                beta,
                self.time_coeffs,
                dirichlet_bcs,
                forcing_zones,
            )

            self.solver = linear_solver_from_input(
                sim, 'solver/rho', SOLVER, PRECONDITIONER, None, KRYLOV_PARAMETERS
            )
            self.slope_limiter = SlopeLimiter(sim, 'rho', self.rho)
            self.slope_limiter.set_phi_old(self.rho_p)

        # Make sure the initial value is included in XDMF results from timestep 0
        self.rho.assign(self.rho_p)

    def get_density(self, k):
        """
        Return the function as defined on timestep t^{n+k}
        """
        if k == 0:
            return self.rho
        elif k == -1:
            return self.rho_p
        elif k == -2:
            return self.rho_pp

    def get_laminar_kinematic_viscosity(self, k):
        """
        It is assumed that the kinematic viscosity is constant
        """
        return Constant(self.nu)

    def get_laminar_dynamic_viscosity(self, k):
        """
        Calculate the blended dynamic viscosity function as a function
        of the (constant) nu and (variable) rho

        Return the function as defined on timestep t^{n+k}
        """
        nu = self.get_laminar_kinematic_viscosity(k)
        rho = self.get_density(k)
        return nu * rho

    def get_density_range(self):
        """
        Return the maximum and minimum densities, rho
        """
        return self.rho_min, self.rho_max

    def get_laminar_kinematic_viscosity_range(self):
        """
        Return the maximum and minimum kinematic viscosities, nu
        """
        return self.nu, self.nu

    def get_laminar_dynamic_viscosity_range(self):
        """
        The minimum and maximum laminar dynamic viscosities, mu.
        """
        return self.nu * self.rho_min, self.nu * self.rho_max

    def update(self, timestep_number, t, dt):
        """
        Update the density field by advecting it for a time dt
        using the given divergence free velocity field
        """
        timer = dolfin.Timer('Ocellaris update rho')
        sim = self.simulation

        if timestep_number != 1:
            # Update the previous values
            self.rho_pp.assign(self.rho_p)
            self.rho_p.assign(self.rho)

        # Check for steady solution every timestep, this can change over time
        force_static = sim.input.get_value(
            'multiphase_solver/force_static', FORCE_STATIC, 'bool'
        )

        if force_static:
            # Keep the existing solution
            self.rho.assign(self.rho_p)

        elif self.use_analytical_solution:
            # Use an analytical density field for testing other parts of Ocellaris
            cpp_code = sim.input.get_value(
                'initial_conditions/rho_p/cpp_code', required_type='string'
            )
            description = 'initial condition for rho_p'
            V = sim.data['Vrho']
            ocellaris_interpolate(sim, cpp_code, description, V, self.rho)

        elif self.use_rk_method:
            # Strong-Stability-Preserving Runge-Kutta DG time integration
            self.rho.assign(self.rho_p)
            self.rho_explicit.assign(self.rho_p)
            self.rk.step(dt)

        else:
            # Compute global bounds
            if self.is_first_timestep:
                lo, hi = self.slope_limiter.set_global_bounds(self.rho)
                if self.slope_limiter.has_global_bounds:
                    sim.log.info(
                        'Setting global bounds [%r, %r] in VariableDensity' % (lo, hi)
                    )

            # Solve the implicit advection equation
            A = self.eq.assemble_lhs()
            b = self.eq.assemble_rhs()
            self.solver.solve(A, self.rho.vector(), b)
            self.slope_limiter.run()
            self.time_coeffs.assign(Constant([3 / 2, -2, 1 / 2]))

        sim.reporting.report_timestep_value('min(rho)', self.rho.vector().min())
        sim.reporting.report_timestep_value('max(rho)', self.rho.vector().max())

        timer.stop()  # Stop timer before hook
        self.simulation.hooks.run_custom_hook('MultiPhaseModelUpdated')
        self.is_first_timestep = False
