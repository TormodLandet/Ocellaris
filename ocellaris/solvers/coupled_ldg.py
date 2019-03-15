# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from dolfin import dx, dS, div, grad, dot, inner, outer, jump, avg, Constant
from ocellaris.utils import ocellaris_error, timeit, linear_solver_from_input
from . import Solver, register_solver, BDM
from ..solver_parts import VelocityBDMProjection
from .coupled import get_global_row_number


# Default values, can be changed in the input file
LU_SOLVER_1CPU = 'default'
LU_SOLVER_NCPU = 'superlu_dist'
LU_PARAMETERS = {}

# Default values, can be changed in the input file
FIX_PRESSURE_DOF = True


@register_solver('CoupledLDG')
class SolverCoupledLDG(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes coupled solver based on the LDG scheme

        :type simulation: ocellaris.Simulation
        """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()

        # First time step timestepping coefficients
        self.set_timestepping_coefficients([1, -1, 0])

        # Solver control parameters
        sim.data['dt'] = Constant(simulation.dt)

        # Create equation
        self.eqs = CoupledEquationsLDG(simulation)

        # Velocity post_processing
        self.velocity_postprocessor = None
        if self.velocity_postprocessing_method == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(sim, sim.data['u'])

        if self.fix_pressure_dof:
            pdof = get_global_row_number(self.subspaces[sim.ndim])
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

        # Give warning if using iterative solver
        if isinstance(self.coupled_solver, dolfin.PETScKrylovSolver):
            sim.log.warning(
                'WARNING: Using a Krylov solver for the coupled NS equations is not a good idea'
            )
        else:
            # Removed in DOLFIN 2018.1
            # self.coupled_solver.set_parameter('same_nonzero_pattern', True)
            pass

        # Deal with pressure null space
        self.fix_pressure_dof = sim.input.get_value(
            'solver/fix_pressure_dof', FIX_PRESSURE_DOF, 'bool'
        )
        # No need for any tricks if the pressure is set via Dirichlet conditions somewhere
        if self.simulation.data['dirichlet_bcs'].get('p', []):
            self.fix_pressure_dof = False

        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        Vp_family = sim.data['Vp'].ufl_element().family()
        if not Vu_family == Vp_family == 'Discontinuous Lagrange':
            ocellaris_error(
                'Must use DG function spaces',
                'LDG solver requires DG spaces for velocity and pressure',
            )

        # Local DG velocity postprocessing
        self.velocity_postprocessing_method = sim.input.get_value(
            'solver/velocity_postprocessing', BDM, 'string'
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
        cd = sim.data['constrained_domain']

        # Create coupled mixed function space and mixed function to hold results
        func_spaces = [Vu] * sim.ndim + [Vp]
        self.subspace_names = ['u%d' % d for d in range(sim.ndim)] + ['p']

        # Create stress tensor space
        P = Vu.ufl_element().degree()
        Vs = dolfin.FunctionSpace(sim.data['mesh'], 'DG', P, constrained_domain=cd)
        for i in range(sim.ndim ** 2):
            stress_name = 'stress_%d' % i
            sim.data[stress_name] = dolfin.Function(Vs)
            func_spaces.append(Vs)
            self.subspace_names.append(stress_name)

        # Create mixed space
        e_mixed = dolfin.MixedElement([fs.ufl_element() for fs in func_spaces])
        Vcoupled = dolfin.FunctionSpace(sim.data['mesh'], e_mixed)
        sim.data['Vcoupled'] = Vcoupled

        # Create function assigner
        Nspace = len(func_spaces)
        self.subspaces = [Vcoupled.sub(i) for i in range(Nspace)]
        sim.data['coupled'] = self.coupled_func = dolfin.Function(Vcoupled)
        self.assigner = dolfin.FunctionAssigner(func_spaces, Vcoupled)

        # Create segregated functions on component and vector form
        u_list, up_list, upp_list, u_conv = [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            u_list.append(u)
            up_list.append(up)
            upp_list.append(upp)
            u_conv.append(uc)
        sim.data['u'] = dolfin.as_vector(u_list)
        sim.data['up'] = dolfin.as_vector(up_list)
        sim.data['upp'] = dolfin.as_vector(upp_list)
        sim.data['u_conv'] = dolfin.as_vector(u_conv)
        sim.data['p'] = dolfin.Function(Vp)

    def set_timestepping_coefficients(self, coeffs):
        """
        Set the time stepping coefficients used for the temporal derivative
        """
        if not 'time_coeffs' in self.simulation.data:
            self.is_first_timestep = True
            self.simulation.data['time_coeffs'] = Constant(coeffs)
            self.simulation.data['time_coeffs_py'] = coeffs
        else:
            self.simulation.data['time_coeffs'].assign(Constant(coeffs))
            self.simulation.data['time_coeffs_py'] = coeffs

    @timeit
    def update_convection(self, t, dt):
        """
        Update terms used to linearise and discretise the convective term
        """
        ndim = self.simulation.ndim
        data = self.simulation.data

        # Update convective velocity field components
        for d in range(ndim):
            uic = data['u_conv%d' % d]
            uip = data['up%d' % d]
            uipp = data['upp%d' % d]

            if self.is_first_timestep:
                uic.assign(uip)
            else:
                # Backwards difference formulation - standard linear extrapolation
                uic.vector().zero()
                uic.vector().axpy(2.0, uip.vector())
                uic.vector().axpy(-1.0, uipp.vector())

        self.is_first_timestep = False

    @timeit
    def postprocess_velocity(self):
        """
        Apply a post-processing operator to the given velocity field
        """
        if self.velocity_postprocessor:
            self.velocity_postprocessor.run()

    @timeit
    def solve_coupled(self):
        """
        Solve the coupled equations
        """
        # Assemble the equation system
        A = self.eqs.assemble_lhs()
        b = self.eqs.assemble_rhs()

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

        # If we remove the null space of the matrix system this will not be the exact same as
        # removing the proper null space of the equation, so we fix this here
        if self.fix_pressure_dof:
            p = self.simulation.data['p']
            dx2 = dolfin.dx(domain=p.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1) * dx2)
            # Perform correction multiple times due to round-of error. The first correction
            # can be i.e 1e14 while the next correction is around unity
            pavg = 1e10
            while abs(pavg) > 1000:
                pavg = dolfin.assemble(p * dx2) / vol
                p.vector()[:] -= pavg

    @timeit
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation
        sim.hooks.simulation_started()
        t = sim.time
        it = sim.timestep

        # Check if there are non-zero values in the upp vectors
        maxabs = 0
        for d in range(sim.ndim):
            this_maxabs = abs(sim.data['upp%d' % d].vector().get_local()).max()
            maxabs = max(maxabs, this_maxabs)
        maxabs = dolfin.MPI.max(dolfin.MPI.comm_world, float(maxabs))
        has_upp_start_values = maxabs > 0

        # Previous-previous values are provided so we can start up with second order time stepping
        if has_upp_start_values:
            sim.log.info('Initial values for upp are found and used')
            self.is_first_timestep = False
            self.set_timestepping_coefficients([3 / 2, -2, 1 / 2])

        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')

            # Check if the simulation is done
            if t + dt > tmax + 1e-6:
                break

            # Advance one time step
            it += 1
            t += dt
            self.simulation.data['dt'].assign(dt)
            self.simulation.hooks.new_timestep(it, t, dt)

            # Extrapolate the convecting velocity to the new time step
            self.update_convection(t, dt)

            # Solve for the new time step
            self.solve_coupled()

            # Postprocess the solution velocity field
            self.postprocess_velocity()

            # Move u -> up, up -> upp and prepare for the next time step
            vel_diff = 0
            for d in range(self.simulation.ndim):
                u_new = self.simulation.data['u%d' % d]
                up = self.simulation.data['up%d' % d]
                upp = self.simulation.data['upp%d' % d]

                if self.is_steady:
                    diff = abs(u_new.vector().get_local() - up.vector().get_local()).max()
                    vel_diff = max(vel_diff, diff)

                upp.assign(up)
                up.assign(u_new)

            # Change time coefficient to second order
            self.set_timestepping_coefficients([3 / 2, -2, 1 / 2])

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


class CoupledEquationsLDG(object):
    def __init__(self, simulation):
        """
        This class assembles the coupled Navier-Stokes equations with LDG discretization

        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
        self.use_grad_p_form = False
        self.use_grad_q_form = True
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
        sigma, tau = [], []
        for d in range(ndim):
            ulist.append(uc[d])
            vlist.append(vc[d])
            indices = list(range(1 + ndim * (d + 1), 1 + ndim * (d + 2)))
            sigma.append(dolfin.as_vector([uc[i] for i in indices]))
            tau.append(dolfin.as_vector([vc[i] for i in indices]))

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

        # Upwind and downwind velocities
        w_nU = (dot(u_conv, n) + abs(dot(u_conv, n))) / 2.0
        w_nD = (dot(u_conv, n) - abs(dot(u_conv, n))) / 2.0

        # LDG penalties
        C11 = Constant(1.0 + sim.ndim)
        switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), n('+'), n('-'))
        C12 = 0.5 * switch

        # Start building the coupled equations
        eq = 0

        # Momentum equations
        for d in range(sim.ndim):
            up = sim.data['up%d' % d]
            upp = sim.data['upp%d' % d]

            # Divergence free criterion
            # ∇⋅u = 0
            u_hat_p = avg(u[d])
            if self.use_grad_q_form:
                eq -= u[d] * q.dx(d) * dx
                eq += u_hat_p * jump(q) * n[d]('+') * dS
            else:
                eq += q * u[d].dx(d) * dx
                eq -= avg(q) * jump(u[d]) * n[d]('+') * dS

            # Time derivative
            # ∂(ρu)/∂t
            eq += rho * (c1 * u[d] + c2 * up + c3 * upp) / dt * v[d] * dx

            # Convection:
            # -w⋅∇(ρu)
            flux_nU = u[d] * w_nU
            flux = jump(flux_nU)
            eq -= u[d] * dot(grad(rho * v[d]), u_conv) * dx
            eq += flux * jump(rho * v[d]) * dS

            # Diffusion:
            # -∇⋅∇u
            u_hat_dS = avg(u[d]) - dot(C12, jump(u[d], n))
            sigma_hat_dS = avg(sigma[d]) - C11 * jump(u[d], n) + C12 * jump(sigma[d], n)
            eq += dot(sigma[d], tau[d]) * dx
            eq += u[d] * div(mu * tau[d]) * dx
            eq -= u_hat_dS * jump(mu * tau[d], n) * dS
            eq += dot(sigma[d], grad(v[d])) * dx
            eq -= dot(sigma_hat_dS, jump(v[d], n)) * dS

            # Pressure
            # ∇p
            if self.use_grad_p_form:
                eq += v[d] * p.dx(d) * dx
                eq -= avg(v[d]) * jump(p) * n[d]('+') * dS
            else:
                eq -= p * v[d].dx(d) * dx
                eq += avg(p) * jump(v[d]) * n[d]('+') * dS

            # Body force (gravity)
            # ρ g
            eq -= rho * g[d] * v[d] * dx

            # Other sources
            for f in sim.data['momentum_sources']:
                eq -= f[d] * v[d] * dx

            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()

                # Divergence free criterion
                if self.use_grad_q_form:
                    eq += q * u_bc * n[d] * dbc.ds()
                else:
                    eq -= q * u[d] * n[d] * dbc.ds()
                    eq += q * u_bc * n[d] * dbc.ds()

                # Convection
                eq += rho * u[d] * w_nU * v[d] * dbc.ds()
                eq += rho * u_bc * w_nD * v[d] * dbc.ds()

                # Diffusion
                u_hat_ds = u_bc
                sigma_hat_ds = sigma[d] - C11 * (u[d] - u_bc) * n
                eq -= u_hat_ds * mu * dot(tau[d], n) * dbc.ds()
                eq -= dot(sigma_hat_ds, n) * v[d] * dbc.ds()

                # Pressure
                if not self.use_grad_p_form:
                    eq += p * v[d] * n[d] * dbc.ds()

            # Neumann boundary
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
            for nbc in neumann_bcs:
                # Divergence free criterion
                if self.use_grad_q_form:
                    eq += q * u[d] * n[d] * nbc.ds()
                else:
                    eq -= q * u[d] * n[d] * nbc.ds()

                # Convection
                eq += rho * u[d] * w_nU * v[d] * nbc.ds()

                # Diffusion
                u_hat_ds = u[d]
                sigma_hat_ds = nbc.func() / mu * n
                eq -= u_hat_ds * mu * dot(tau[d], n) * nbc.ds()
                eq -= dot(sigma_hat_ds, n) * v[d] * nbc.ds()

                # Pressure
                if not self.use_grad_p_form:
                    eq += p * v[d] * n[d] * nbc.ds()

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


class CoupledEquationsLDG2(object):
    def __init__(self, simulation):
        """
        This class assembles the coupled Navier-Stokes equations with LDG discretization

        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
        self.use_grad_p_form = False
        self.use_grad_q_form = True
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
        sigmas, taus = [], []
        for d in range(sim.ndim):
            ulist.append(uc[d])
            vlist.append(vc[d])
            indices = list(range(1 + sim.ndim * (d + 1), 1 + sim.ndim * (d + 2)))
            sigmas.append([uc[i] for i in indices])
            taus.append([vc[i] for i in indices])

        u = dolfin.as_vector(ulist)
        v = dolfin.as_vector(vlist)
        p = uc[sim.ndim]
        q = vc[sim.ndim]
        sigma = dolfin.as_tensor(sigmas)
        tau = dolfin.as_tensor(taus)

        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        n = dolfin.FacetNormal(mesh)
        h = dolfin.FacetArea(mesh)

        # Fluid properties
        rho = mpm.get_density(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)

        # Upwind and downwind velocities
        w_nU = (dot(u_conv, n) + abs(dot(u_conv, n))) / 2.0
        w_nD = (dot(u_conv, n) - abs(dot(u_conv, n))) / 2.0
        u_uw_s = dolfin.conditional(dolfin.gt(dot(u_conv, n), 0.0), 1.0, 0.0)('+')
        u_uw = u_uw_s * u('+') + (1 - u_uw_s) * u('-')

        # LDG penalties
        # kappa_0 = Constant(4.0)
        # kappa = mu*kappa_0/h
        C11 = avg(mu / h)
        D11 = avg(h / mu)
        C12 = 0.2 * n('+')
        D12 = 0.2 * n('+')

        def ojump(v, n):
            return outer(v, n)('+') + outer(v, n)('-')

        # Interior facet fluxes
        # u_hat_dS = avg(u)
        # sigma_hat_dS = avg(sigma) - avg(kappa)*ojump(u, n)
        # p_hat_dS = avg(p)
        u_hat_s_dS = avg(u) + dot(ojump(u, n), C12)
        u_hat_p_dS = avg(u) + D11 * jump(p, n) + D12 * jump(u, n)
        sigma_hat_dS = avg(sigma) - C11 * ojump(u, n) - outer(jump(sigma, n), C12)
        p_hat_dS = avg(p) - dot(D12, jump(p, n))

        # Time derivative
        up = sim.data['up']
        upp = sim.data['upp']
        eq = rho * dot(c1 * u + c2 * up + c3 * upp, v) / dt * dx

        # LDG equation 1
        eq += inner(sigma, tau) * dx
        eq += dot(u, div(mu * tau)) * dx
        eq -= dot(u_hat_s_dS, jump(mu * tau, n)) * dS

        # LDG equation 2
        eq += (inner(sigma, grad(v)) - p * div(v)) * dx
        eq -= (inner(sigma_hat_dS, ojump(v, n)) - p_hat_dS * jump(v, n)) * dS
        eq -= dot(u, div(outer(v, rho * u_conv))) * dx
        eq += rho('+') * dot(u_conv('+'), n('+')) * dot(u_uw, v('+')) * dS
        eq += rho('-') * dot(u_conv('-'), n('-')) * dot(u_uw, v('-')) * dS
        momentum_sources = sim.data['momentum_sources'] + [rho * g]
        eq -= dot(sum(momentum_sources), v) * dx

        # LDG equation 3
        eq -= dot(u, grad(q)) * dx
        eq += dot(u_hat_p_dS, jump(q, n)) * dS

        # Dirichlet boundary
        dirichlet_bcs = get_collected_velocity_bcs(sim, 'dirichlet_bcs')
        for ds, u_bc in dirichlet_bcs.items():
            # sigma_hat_ds = sigma - kappa*outer(u, n)
            sigma_hat_ds = sigma - C11 * outer(u - u_bc, n)
            u_hat_ds = u_bc
            p_hat_ds = p

            # LDG equation 1
            eq -= dot(u_hat_ds, dot(mu * tau, n)) * ds

            # LDG equation 2
            eq -= (inner(sigma_hat_ds, outer(v, n)) - p_hat_ds * dot(v, n)) * ds
            eq += rho * w_nU * dot(u, v) * ds
            eq += rho * w_nD * dot(u_bc, v) * ds

            # LDG equation 3
            eq += dot(u_hat_ds, q * n) * ds

        # Neumann boundary
        neumann_bcs = get_collected_velocity_bcs(sim, 'neumann_bcs')
        assert not neumann_bcs
        for ds, du_bc in neumann_bcs.items():
            # Divergence free criterion
            if self.use_grad_q_form:
                eq += q * dot(u, n) * ds
            else:
                eq -= q * dot(u, n) * ds

            # Convection
            eq += rho * w_nU * dot(u, v) * ds

            # Diffusion
            u_hat_ds = u
            sigma_hat_ds = outer(du_bc, n) / mu
            eq -= dot(u_hat_ds, dot(mu * tau, n)) * ds
            eq -= inner(sigma_hat_ds, outer(v, n)) * ds

            # Pressure
            if not self.use_grad_p_form:
                eq += p * dot(v, n) * ds

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


def get_collected_velocity_bcs(simulation, name):
    """
    When mixed Dirichlet/Neumann BCs on the same facet (for different velocity
    components) is not supported it can be convenient to get all BCs collected
    by boundary region. This function returns a dictionary for Dirichlet or
    Neumann BCs where the "ds" meassure is the key and a vector of boundary
    values for each velocity component is the value

    The "name" parameter must be 'dirichlet_bcs' or 'neumann_bcs'
    """
    # Collect BCs
    bc_dict = {}
    for d in range(simulation.ndim):
        for bc in simulation.data[name].get('u%d' % d, []):
            if d == 0:
                bc_dict[bc.ds()] = [bc.func()]
            else:
                bc_dict[bc.ds()].append(bc.func())

    # Verify that all components are present and convert to vector
    for ds, u_bc in bc_dict.items():
        assert len(u_bc) == simulation.ndim
        bc_dict[ds] = dolfin.as_vector(u_bc)

    simulation.log.debug('    Found %d %s boundary regions' % (len(bc_dict), name))
    return bc_dict
