# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ufl.classes import Zero
from dolfin import dot, div, grad, jump


class AdvectionEquation(object):
    def __init__(
        self,
        simulation,
        Vc,
        cp,
        cpp,
        u_conv,
        beta,
        time_coeffs,
        dirichlet_bcs,
        forcing_zones=(),
        dt=None,
    ):
        """
        This class assembles the advection equation for a scalar function c
        """
        self.simulation = simulation
        self.Vc = Vc
        self.cp = cp
        self.cpp = cpp
        self.u_conv = u_conv
        self.beta = beta
        self.time_coeffs = time_coeffs
        self.dirichlet_bcs = dirichlet_bcs
        self.forcing_zones = forcing_zones
        self.dt = dt if dt is not None else simulation.data['dt']

        # Discontinuous or continuous elements
        Vc_family = Vc.ufl_element().family()
        self.colour_is_discontinuous = Vc_family == 'Discontinuous Lagrange'

        if isinstance(u_conv[0], (int, float, dolfin.Constant, Zero)):
            self.velocity_is_trace = False
        else:
            try:
                Vu_family = u_conv[0].function_space().ufl_element().family()
            except Exception:
                Vu_family = u_conv.function_space().ufl_element().family()
            self.velocity_is_trace = Vu_family == 'Discontinuous Lagrange Trace'

        if self.velocity_is_trace:
            assert self.colour_is_discontinuous
            assert Vc.ufl_element().degree() == 0

        # Create UFL forms
        self.define_advection_equation()

    def define_advection_equation(self):
        """
        Setup the advection equation for the colour function

        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        dS, dx = dolfin.dS(mesh), dolfin.dx(mesh)

        # Trial and test functions
        Vc = self.Vc
        c = dolfin.TrialFunction(Vc)
        d = dolfin.TestFunction(Vc)

        c1, c2, c3 = self.time_coeffs
        dt = self.dt
        u_conv = self.u_conv

        if not self.colour_is_discontinuous:
            # Continous Galerkin implementation of the advection equation
            # FIXME: add stabilization
            eq = (c1 * c + c2 * self.cp + c3 * self.cpp) / dt * d * dx + div(c * u_conv) * d * dx

        elif self.velocity_is_trace:
            # Upstream and downstream normal velocities
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n))) / 2
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n))) / 2

            if self.beta is not None:
                # Define the blended flux
                # The blending factor beta is not DG, so beta('+') == beta('-')
                b = self.beta('+')
                flux = (1 - b) * jump(c * w_nU) + b * jump(c * w_nD)
            else:
                flux = jump(c * w_nU)

            # Discontinuous Galerkin implementation of the advection equation
            eq = (c1 * c + c2 * self.cp + c3 * self.cpp) / dt * d * dx + flux * jump(d) * dS

            # On each facet either w_nD or w_nU will be 0, the other is multiplied
            # with the appropriate flux, either the value c going out of the domain
            # or the Dirichlet value coming into the domain
            for dbc in self.dirichlet_bcs:
                eq += w_nD * dbc.func() * d * dbc.ds()
                eq += w_nU * c * d * dbc.ds()

        elif self.beta is not None:
            # Upstream and downstream normal velocities
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n))) / 2
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n))) / 2

            if self.beta is not None:
                # Define the blended flux
                # The blending factor beta is not DG, so beta('+') == beta('-')
                b = self.beta('+')
                flux = (1 - b) * jump(c * w_nU) + b * jump(c * w_nD)
            else:
                flux = jump(c * w_nU)

            # Discontinuous Galerkin implementation of the advection equation
            eq = (
                (c1 * c + c2 * self.cp + c3 * self.cpp) / dt * d * dx
                - dot(c * u_conv, grad(d)) * dx
                + flux * jump(d) * dS
            )

            # Enforce Dirichlet BCs weakly
            for dbc in self.dirichlet_bcs:
                eq += w_nD * dbc.func() * d * dbc.ds()
                eq += w_nU * c * d * dbc.ds()

        else:
            # Downstream normal velocities
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n))) / 2

            # Discontinuous Galerkin implementation of the advection equation
            eq = (c1 * c + c2 * self.cp + c3 * self.cpp) / dt * d * dx

            # Convection integrated by parts two times to bring back the original
            # div form (this means we must subtract and add all fluxes)
            eq += div(c * u_conv) * d * dx

            # Replace downwind flux with upwind flux on downwind internal facets
            eq -= jump(w_nD * d) * jump(c) * dS

            # Replace downwind flux with upwind BC flux on downwind external facets
            for dbc in self.dirichlet_bcs:
                # Subtract the "normal" downwind flux
                eq -= w_nD * c * d * dbc.ds()
                # Add the boundary value upwind flux
                eq += w_nD * dbc.func() * d * dbc.ds()

        # Penalty forcing zones
        for fz in self.forcing_zones:
            eq += fz.penalty * fz.beta * (c - fz.target) * d * dx

        a, L = dolfin.system(eq)
        self.form_lhs = dolfin.Form(a)
        self.form_rhs = dolfin.Form(L)
        self.tensor_lhs = None
        self.tensor_rhs = None

    def assemble_lhs(self):
        if self.tensor_lhs is None:
            lhs = dolfin.assemble(self.form_lhs)
            self.tensor_lhs = dolfin.as_backend_type(lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            rhs = dolfin.assemble(self.form_rhs)
            self.tensor_rhs = dolfin.as_backend_type(rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs
