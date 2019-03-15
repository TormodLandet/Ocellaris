# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from dolfin import dx, div, grad, dot


class CoupledEquationsCG(object):
    use_strong_bcs = True

    def __init__(
        self,
        simulation,
        flux_type,
        use_stress_divergence_form,
        use_grad_p_form,
        use_grad_q_form,
        use_lagrange_multiplicator,
        pressure_continuity_factor,
        velocity_continuity_factor_D12,
        include_hydrostatic_pressure,
        incompressibility_flux_type,
    ):
        """
        Weak form of the Navier-Stokes eq. on coupled form with continuous elements

        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.use_grad_q_form = use_grad_q_form
        self.flux_type = flux_type
        self.use_lagrange_multiplicator = use_lagrange_multiplicator
        self.pressure_continuity_factor = pressure_continuity_factor
        self.velocity_continuity_factor_D12 = velocity_continuity_factor_D12
        self.include_hydrostatic_pressure = include_hydrostatic_pressure
        self.incompressibility_flux_type = incompressibility_flux_type

        assert self.incompressibility_flux_type in ('central', 'upwind')

        # Create UFL forms
        self.define_coupled_equation()

    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation

        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mpm = sim.multi_phase_model
        mesh = sim.data['mesh']
        Vcoupled = sim.data['Vcoupled']
        u_conv = sim.data['u_conv']

        # Unpack the coupled trial and test functions
        uc = dolfin.TrialFunction(Vcoupled)
        vc = dolfin.TestFunction(Vcoupled)
        ulist = []
        vlist = []
        ndim = self.simulation.ndim
        for d in range(ndim):
            ulist.append(uc[d])
            vlist.append(vc[d])

        u = dolfin.as_vector(ulist)
        v = dolfin.as_vector(vlist)
        p = uc[ndim]
        q = vc[ndim]

        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        n = dolfin.FacetNormal(mesh)

        # Fluid properties
        rho = mpm.get_density(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)

        # Hydrostatic pressure correction
        if self.include_hydrostatic_pressure:
            p += sim.data['p_hydrostatic']

        # Start building the coupled equations
        eq = 0

        # ALE mesh velocities
        if sim.mesh_morpher.active:
            u_mesh = sim.data['u_mesh']

            # Either modify the convective velocity or just include the mesh
            # velocity on cell integral form. Only activate one of these lines
            # u_conv -= u_mesh
            eq -= dot(div(rho * dolfin.outer(u, u_mesh)), v) * dx

            # Divergence of u should balance expansion/contraction of the cell K
            # ∇⋅u = -∂x/∂t       (See below for definition of the ∇⋅u term)
            cvol_new = dolfin.CellVolume(mesh)
            cvol_old = sim.data['cvolp']
            eq += (cvol_new - cvol_old) / dt * q * dx

        # Lagrange multiplicator to remove the pressure null space
        # ∫ p dx = 0
        if self.use_lagrange_multiplicator:
            lm_trial = uc[ndim + 1]
            lm_test = vc[ndim + 1]
            eq = (p * lm_test + q * lm_trial) * dx

        # Momentum equations
        for d in range(sim.ndim):
            up = sim.data['up%d' % d]
            upp = sim.data['upp%d' % d]

            # Divergence free criterion
            # ∇⋅u = 0
            eq += u[d].dx(d) * q * dx

            # Time derivative
            # ∂u/∂t
            eq += rho * (c1 * u[d] + c2 * up + c3 * upp) / dt * v[d] * dx

            # Convection
            # ∇⋅(ρ u ⊗ u_conv)
            eq += rho * dot(u_conv, grad(u[d])) * v[d] * dx

            # Diffusion
            # -∇⋅μ(∇u)
            eq += mu * dot(grad(u[d]), grad(v[d])) * dx

            # -∇⋅μ(∇u)^T
            if self.use_stress_divergence_form:
                eq += mu * dot(u.dx(d), grad(v[d])) * dx

            # Pressure
            # ∇p
            eq -= v[d].dx(d) * p * dx

            # Body force (gravity)
            # ρ g
            eq -= rho * g[d] * v[d] * dx

            # Other sources
            for f in sim.data['momentum_sources']:
                eq -= f[d] * v[d] * dx

            # Neumann boundary conditions
            neumann_bcs_pressure = sim.data['neumann_bcs'].get('p', [])
            for nbc in neumann_bcs_pressure:
                eq += p * v[d] * n[d] * nbc.ds()

            # Outlet boundary
            for obc in sim.data['outlet_bcs']:
                # Diffusion
                mu_dudn = p * n[d]
                eq -= mu_dudn * v[d] * obc.ds()

                # Pressure
                p_ = mu * dot(dot(grad(u), n), n)
                eq += p_ * v[d] * n[d] * obc.ds()

        a, L = dolfin.system(eq)
        self.form_lhs = a
        self.form_rhs = L
        self.tensor_lhs = None
        self.tensor_rhs = None

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
