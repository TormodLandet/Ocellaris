# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from dolfin import FiniteElement, VectorElement, MixedElement, FunctionSpace, VectorFunctionSpace
from dolfin import FacetNormal, TrialFunction, TestFunction, TestFunctions, Function
from dolfin import dot, as_vector, dx, dS, ds, LocalSolver


class VelocityBDMProjection:
    def __init__(
        self,
        simulation,
        w,
        incompressibility_flux_type='central',
        D12=None,
        degree=None,
        use_bcs=True,
        use_nedelec=True,
    ):
        """
        Implement equation 4a and 4b in "Two new techniques for generating exactly
        incompressible approximate velocities" by Bernardo Cocburn (2009)

        For each element K in the mesh:

            <u⋅n, φ> = <û⋅n, φ>  ∀ ϕ ∈ P_{k}(F) for any face F ∈ ∂K
            (u, ϕ) = (w, ϕ)      ∀ φ ∈ P_{k-2}(K)^2
            (u, ϕ) = (w, ϕ)      ∀ φ ∈ {ϕ ∈ P_{k}(K)^2 : ∇⋅ϕ = 0 in K, ϕ⋅n = 0 on ∂K}

        Here w is the input velocity function in DG2 space and û is the flux at
        each face. P_{x} is the space of polynomials of order k

        The flux type can be 'central' or 'upwind'
        """
        self.simulation = simulation
        simulation.log.info('    Setting up velocity BDM projection')

        V = w[0].function_space()
        ue = V.ufl_element()
        gdim = w.ufl_shape[0]
        if degree is None:
            pdeg = ue.degree()
            Vout = V
        else:
            pdeg = degree
            Vout = FunctionSpace(V.mesh(), 'DG', degree)
        pg = (pdeg, gdim)

        assert ue.family() == 'Discontinuous Lagrange'
        assert incompressibility_flux_type in ('central', 'upwind')

        if use_nedelec and pdeg > 1:
            a, L, V = self._setup_projection_nedelec(
                w, incompressibility_flux_type, D12, use_bcs, pdeg, gdim
            )
        elif gdim == 2 and pdeg == 1:
            a, L, V = self._setup_dg1_projection_2D(w, incompressibility_flux_type, D12, use_bcs)
        elif gdim == 2 and pdeg == 2:
            a, L, V = self._setup_dg2_projection_2D(w, incompressibility_flux_type, D12, use_bcs)
        else:
            raise NotImplementedError(
                'VelocityBDMProjection does not support ' 'degree %d and dimension %d' % pg
            )

        # Pre-factorize matrices and store for usage in projection
        self.local_solver = LocalSolver(a, L)
        self.local_solver.factorize()
        self.temp_function = Function(V)
        self.w = w

        # Create function assigners
        self.assigners = []
        for i in range(gdim):
            self.assigners.append(dolfin.FunctionAssigner(Vout, V.sub(i)))

    def _setup_dg1_projection_2D(self, w, incompressibility_flux_type, D12, use_bcs):
        """
        Implement the projection where the result is BDM embeded in a DG1 function
        """
        sim = self.simulation
        k = 1
        gdim = 2
        mesh = w[0].function_space().mesh()
        V = VectorFunctionSpace(mesh, 'DG', k)
        W = FunctionSpace(mesh, 'DGT', k)
        n = FacetNormal(mesh)

        v1 = TestFunction(W)
        u = TrialFunction(V)

        # The same fluxes that are used in the incompressibility equation
        if incompressibility_flux_type == 'central':
            u_hat_dS = dolfin.avg(w)
        elif incompressibility_flux_type == 'upwind':
            w_nU = (dot(w, n) + abs(dot(w, n))) / 2.0
            switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), 1.0, 0.0)
            u_hat_dS = switch * w('+') + (1 - switch) * w('-')

        if D12 is not None:
            u_hat_dS += dolfin.Constant([D12, D12]) * dolfin.jump(w, n)

        # Equation 1 - flux through the sides
        a = L = 0
        for R in '+-':
            a += dot(u(R), n(R)) * v1(R) * dS
            L += dot(u_hat_dS, n(R)) * v1(R) * dS

        # Eq. 1 cont. - flux through external boundaries
        a += dot(u, n) * v1 * ds
        if use_bcs:
            for d in range(gdim):
                dirichlet_bcs = sim.data['dirichlet_bcs']['u%d' % d]
                neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
                robin_bcs = sim.data['robin_bcs'].get('u%d' % d, [])
                outlet_bcs = sim.data['outlet_bcs']

                for dbc in dirichlet_bcs:
                    u_bc = dbc.func()
                    L += u_bc * n[d] * v1 * dbc.ds()

                for nbc in neumann_bcs + robin_bcs + outlet_bcs:
                    if nbc.enforce_zero_flux:
                        pass  # L += 0
                    else:
                        L += w[d] * n[d] * v1 * nbc.ds()

            for sbc in sim.data['slip_bcs'].get('u', []):
                pass  # L += 0
        else:
            L += dot(w, n) * v1 * ds

        # Equation 2 - internal shape   :   empty for DG1
        # Equation 3 - BDM Phi          :   empty for DG1

        return a, L, V

    def _setup_dg2_projection_2D(self, w, incompressibility_flux_type, D12, use_bcs):
        """
        Implement the projection where the result is BDM embeded in a DG2 function
        """
        sim = self.simulation
        k = 2
        gdim = 2
        mesh = w[0].function_space().mesh()
        V = VectorFunctionSpace(mesh, 'DG', k)
        n = FacetNormal(mesh)

        # The mixed function space of the projection test functions
        e1 = FiniteElement('DGT', mesh.ufl_cell(), k)
        e2 = VectorElement('DG', mesh.ufl_cell(), k - 2)
        e3 = FiniteElement('Bubble', mesh.ufl_cell(), 3)
        em = MixedElement([e1, e2, e3])
        W = FunctionSpace(mesh, em)
        v1, v2, v3b = TestFunctions(W)
        u = TrialFunction(V)

        # The same fluxes that are used in the incompressibility equation
        if incompressibility_flux_type == 'central':
            u_hat_dS = dolfin.avg(w)
        elif incompressibility_flux_type == 'upwind':
            w_nU = (dot(w, n) + abs(dot(w, n))) / 2.0
            switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), 1.0, 0.0)
            u_hat_dS = switch * w('+') + (1 - switch) * w('-')

        if D12 is not None:
            u_hat_dS += dolfin.Constant([D12, D12]) * dolfin.jump(w, n)

        # Equation 1 - flux through the sides
        a = L = 0
        for R in '+-':
            a += dot(u(R), n(R)) * v1(R) * dS
            L += dot(u_hat_dS, n(R)) * v1(R) * dS

        # Eq. 1 cont. - flux through external boundaries
        a += dot(u, n) * v1 * ds
        if use_bcs:
            for d in range(gdim):
                dirichlet_bcs = sim.data['dirichlet_bcs']['u%d' % d]
                neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
                robin_bcs = sim.data['robin_bcs'].get('u%d' % d, [])
                outlet_bcs = sim.data['outlet_bcs']

                for dbc in dirichlet_bcs:
                    u_bc = dbc.func()
                    L += u_bc * n[d] * v1 * dbc.ds()

                for nbc in neumann_bcs + robin_bcs + outlet_bcs:
                    if nbc.enforce_zero_flux:
                        pass  # L += 0
                    else:
                        L += w[d] * n[d] * v1 * nbc.ds()

            for sbc in sim.data['slip_bcs'].get('u', []):
                pass  # L += 0
        else:
            L += dot(w, n) * v1 * ds

        # Equation 2 - internal shape
        a += dot(u, v2) * dx
        L += dot(w, v2) * dx

        # Equation 3 - BDM Phi
        v3 = as_vector([v3b.dx(1), -v3b.dx(0)])  # Curl of [0, 0, v3b]
        a += dot(u, v3) * dx
        L += dot(w, v3) * dx

        return a, L, V

    def _setup_projection_nedelec(self, w, incompressibility_flux_type, D12, use_bcs, pdeg, gdim):
        """
        Implement the BDM-like projection using Nedelec elements in the test function
        """
        sim = self.simulation
        k = pdeg
        mesh = w[0].function_space().mesh()
        V = VectorFunctionSpace(mesh, 'DG', k)
        n = FacetNormal(mesh)

        # The mixed function space of the projection test functions
        e1 = FiniteElement('DGT', mesh.ufl_cell(), k)
        e2 = FiniteElement('N1curl', mesh.ufl_cell(), k - 1)
        em = MixedElement([e1, e2])
        W = FunctionSpace(mesh, em)
        v1, v2 = TestFunctions(W)
        u = TrialFunction(V)

        # The same fluxes that are used in the incompressibility equation
        if incompressibility_flux_type == 'central':
            u_hat_dS = dolfin.avg(w)
        elif incompressibility_flux_type == 'upwind':
            w_nU = (dot(w, n) + abs(dot(w, n))) / 2.0
            switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), 1.0, 0.0)
            u_hat_dS = switch * w('+') + (1 - switch) * w('-')

        if D12 is not None:
            u_hat_dS += dolfin.Constant([D12] * gdim) * dolfin.jump(w, n)

        # Equation 1 - flux through the sides
        a = L = 0
        for R in '+-':
            a += dot(u(R), n(R)) * v1(R) * dS
            L += dot(u_hat_dS, n(R)) * v1(R) * dS

        # Eq. 1 cont. - flux through external boundaries
        a += dot(u, n) * v1 * ds
        if use_bcs:
            for d in range(gdim):
                dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
                neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
                robin_bcs = sim.data['robin_bcs'].get('u%d' % d, [])
                outlet_bcs = sim.data['outlet_bcs']

                for dbc in dirichlet_bcs:
                    u_bc = dbc.func()
                    L += u_bc * n[d] * v1 * dbc.ds()

                for nbc in neumann_bcs + robin_bcs + outlet_bcs:
                    if nbc.enforce_zero_flux:
                        pass  # L += 0
                    else:
                        L += w[d] * n[d] * v1 * nbc.ds()

            for sbc in sim.data['slip_bcs'].get('u', []):
                pass  # L += 0
        else:
            L += dot(w, n) * v1 * ds

        # Equation 2 - internal shape using 'Nedelec 1st kind H(curl)' elements
        a += dot(u, v2) * dx
        L += dot(w, v2) * dx

        return a, L, V

    def run(self, w=None):
        """
        Perform the projection based on the current state of the Function w
        """
        # Find the projected velocity
        self.local_solver.solve_local_rhs(self.temp_function)

        # Assign to w
        w = self.w if w is None else w
        U = self.temp_function.split()
        for i, a in enumerate(self.assigners):
            a.assign(w[i], U[i])
