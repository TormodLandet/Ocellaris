# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ufl.constantvalue import Zero
from dolfin import dot, grad, avg, jump, dx, dS, Constant
from ..solver_parts import define_penalty
from .coupled_equations import define_dg_equations


class BaseEquation(object):
    # Will be shadowed by object properties after first assemble
    tensor_lhs = None
    tensor_rhs = None

    def assemble_lhs(self):
        if self.tensor_lhs is None:
            self.tensor_lhs = dolfin.assemble(self.form_lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            self.tensor_rhs = dolfin.assemble(self.form_rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs


class MomentumPredictionEquation(BaseEquation):
    def __init__(
        self,
        simulation,
        use_stress_divergence_form,
        use_grad_p_form,
        include_hydrostatic_pressure,
    ):
        """
        This class assembles the standard Navier-Stokes momentum equation
        with an explicit pressure. It uses the coupled form, with ufl Zero
        as pressure test function to only get the (u, v) weak form matrix
        """
        self.simulation = simulation
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.include_hydrostatic_pressure = include_hydrostatic_pressure
        self.define_momentum_equation()

    def define_momentum_equation(self):
        """
        Setup the momentum equation weak form
        """
        sim = self.simulation
        Vuvw = sim.data['uvw_star'].function_space()
        tests = dolfin.TestFunctions(Vuvw)
        trials = dolfin.TrialFunctions(Vuvw)

        # Split into components
        v = dolfin.as_vector(tests[:])
        u = dolfin.as_vector(trials[:])

        # The pressure is explicit p* and q is zero (on a domain, to avoid warnings)
        p = sim.data['p']

        class MyZero(Zero):
            def ufl_domains(self):
                return p.ufl_domains()

        q = MyZero()

        lm_trial = lm_test = None

        # Define the momentum equation weak form
        eq = define_dg_equations(
            u,
            v,
            p,
            q,
            lm_trial,
            lm_test,
            self.simulation,
            include_hydrostatic_pressure=self.include_hydrostatic_pressure,
            incompressibility_flux_type='central',  # Only used with q
            use_grad_q_form=False,  # Only used with q
            use_grad_p_form=self.use_grad_p_form,
            use_stress_divergence_form=self.use_stress_divergence_form,
        )
        self.form_lhs, self.form_rhs = dolfin.system(eq)


class PressureCorrectionEquation(BaseEquation):
    def __init__(
        self, simulation, use_lagrange_multiplicator, incompressibility_flux_type
    ):
        """
        This class assembles the pressure Poisson equation, both CG and DG
        """
        self.simulation = simulation
        self.use_lagrange_multiplicator = use_lagrange_multiplicator
        self.incompressibility_flux_type = incompressibility_flux_type
        self.define_pressure_equation()

    def calculate_penalties(self):
        """
        Calculate SIPG penalty
        """
        mesh = self.simulation.data['mesh']
        P = self.simulation.data['Vp'].ufl_element().degree()
        rho_min, rho_max = self.simulation.multi_phase_model.get_density_range()
        k_min = 1.0 / rho_max
        k_max = 1.0 / rho_min
        penalty_dS = define_penalty(mesh, P, k_min, k_max, boost_factor=3, exponent=1.0)
        penalty_ds = penalty_dS * 2
        self.simulation.log.info(
            '    DG SIP penalty pressure:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds)
        )

        return Constant(penalty_dS), Constant(penalty_ds)

    def define_pressure_equation(self):
        """
        Setup the pressure Poisson equation

        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        Vp = sim.data['Vp']
        p_star = sim.data['p']
        u_star = sim.data['u']

        # Trial and test functions
        p = dolfin.TrialFunction(Vp)
        q = dolfin.TestFunction(Vp)

        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)

        # Fluid properties
        mpm = sim.multi_phase_model
        mu = mpm.get_laminar_dynamic_viscosity(0)
        rho = sim.data['rho']

        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        assert not self.use_lagrange_multiplicator, 'NOT IMPLEMENTED YET'

        # Weak form of the Poisson eq. with discontinuous elements
        # -∇⋅∇p = - γ_1/Δt ρ ∇⋅u^*
        K = 1.0 / rho
        a = K * dot(grad(p), grad(q)) * dx
        L = K * dot(grad(p_star), grad(q)) * dx

        # RHS, ∇⋅u^*
        if self.incompressibility_flux_type == 'central':
            u_flux = avg(u_star)
        elif self.incompressibility_flux_type == 'upwind':
            switch = dolfin.conditional(
                dolfin.gt(abs(dot(u_star, n))('+'), 0.0), 1.0, 0.0
            )
            u_flux = switch * u_star('+') + (1 - switch) * u_star('-')
        L += c1 / dt * dot(u_star, grad(q)) * dx
        L -= c1 / dt * dot(u_flux, n('+')) * jump(q) * dS

        # Symmetric Interior Penalty method for -∇⋅∇p
        a -= dot(n('+'), avg(K * grad(p))) * jump(q) * dS
        a -= dot(n('+'), avg(K * grad(q))) * jump(p) * dS

        # Symmetric Interior Penalty method for -∇⋅∇p^*
        L -= dot(n('+'), avg(K * grad(p_star))) * jump(q) * dS
        L -= dot(n('+'), avg(K * grad(q))) * jump(p_star) * dS

        # Weak continuity
        penalty_dS, penalty_ds = self.calculate_penalties()

        # Symmetric Interior Penalty coercivity term
        a += penalty_dS * jump(p) * jump(q) * dS
        # L += penalty_dS*jump(p_star)*jump(q)*dS

        # Collect Dirichlet and outlet boundary values
        dirichlet_vals_and_ds = []
        for dbc in sim.data['dirichlet_bcs'].get('p', []):
            dirichlet_vals_and_ds.append((dbc.func(), dbc.ds()))
        for obc in sim.data['outlet_bcs']:
            p_ = mu * dot(dot(grad(u_star), n), n)
            dirichlet_vals_and_ds.append((p_, obc.ds()))

        # Apply Dirichlet boundary conditions
        for p_bc, dds in dirichlet_vals_and_ds:
            # SIPG for -∇⋅∇p
            a -= dot(n, K * grad(p)) * q * dds
            a -= dot(n, K * grad(q)) * p * dds
            L -= dot(n, K * grad(q)) * p_bc * dds

            # SIPG for -∇⋅∇p^*
            L -= dot(n, K * grad(p_star)) * q * dds
            L -= dot(n, K * grad(q)) * p_star * dds

            # Weak Dirichlet
            a += penalty_ds * p * q * dds
            L += penalty_ds * p_bc * q * dds

            # Weak Dirichlet for p^*
            # L += penalty_ds*p_star*q*dds
            # L -= penalty_ds*p_bc*q*dds

        # Neumann boundary conditions
        neumann_bcs = sim.data['neumann_bcs'].get('p', [])
        for nbc in neumann_bcs:
            # Neumann boundary conditions on p and p_star cancel
            # L += (nbc.func() - dot(n, grad(p_star)))*q*nbc.ds()
            pass

        # Use boundary conditions for the velocity for the
        # term from integration by parts of div(u_star)
        for d in range(sim.ndim):
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                L -= c1 / dt * u_bc * n[d] * q * dbc.ds()
            for nbc in neumann_bcs:
                L -= c1 / dt * u_star[d] * n[d] * q * nbc.ds()

        # ALE mesh velocities
        if sim.mesh_morpher.active:
            cvol_new = dolfin.CellVolume(mesh)
            cvol_old = sim.data['cvolp']

            # Divergence of u should balance expansion/contraction of the cell K
            # ∇⋅u = -∂K/∂t       (See below for definition of the ∇⋅u term)
            L -= (cvol_new - cvol_old) / dt * q * dx

        self.form_lhs = a
        self.form_rhs = L


class VelocityUpdateEquation(BaseEquation):
    def __init__(self, simulation, component):
        """
        Define the velocity update equation for velocity component d.
        """
        self.simulation = simulation
        self.component = component
        self.define_update_equation()

    def define_update_equation(self):
        sim = self.simulation
        rho = sim.data['rho']
        c1 = sim.data['time_coeffs'][0]
        dt = sim.data['dt']

        Vu = sim.data['Vu']
        us = sim.data['u%d' % self.component]
        p_hat = sim.data['p_hat']
        u = dolfin.TrialFunction(Vu)
        v = dolfin.TestFunction(Vu)

        self.form_lhs = u * v * dx
        self.form_rhs = (
            us * v * dx - dt / (c1 * rho) * p_hat.dx(self.component) * v * dx
        )


EQUATION_SUBTYPES = {
    'Default': (
        MomentumPredictionEquation,
        PressureCorrectionEquation,
        VelocityUpdateEquation,
    )
}
