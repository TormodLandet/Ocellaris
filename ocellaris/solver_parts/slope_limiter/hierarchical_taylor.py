# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin as df
from ocellaris.utils import verify_key, OcellarisError
from ocellaris.utils import lagrange_to_taylor, taylor_to_lagrange, get_local, set_local
from . import register_slope_limiter, SlopeLimiterBase
from .limiter_cpp_utils import SlopeLimiterInput


@register_slope_limiter('HierarchicalTaylor')
class HierarchicalTaylorSlopeLimiter(SlopeLimiterBase):
    description = 'Uses a Taylor DG decomposition to limit derivatives at the vertices'

    def __init__(self, simulation, phi_name, phi, skip_cells, boundary_conditions, limiter_input):
        """
        Limit the slope of the given scalar to obtain boundedness
        """
        # Verify input
        V = phi.function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        tdim = mesh.topology().dim()
        gdim = mesh.geometry().dim()
        loc = 'HierarchalTaylor slope limiter'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (0, 1, 2), loc)
        verify_key('function shape', phi.ufl_shape, [()], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2, 3], loc)
        assert gdim == tdim, (
            'HierarchalTaylor slope limiter requires that '
            'topological and geometrical dimensions are identical'
        )

        # Read input
        use_cpp = limiter_input.get_value('use_cpp', True, 'bool')
        enforce_bounds = limiter_input.get_value('enforce_bounds', False, 'bool')
        enforce_bcs = limiter_input.get_value('enforce_bcs', True, 'bool')
        use_weak_bcs = limiter_input.get_value('use_weak_bcs', True, 'bool')
        trust_robin_dval = limiter_input.get_value('trust_robin_dval', True, 'bool')
        simulation.log.info('        Enforce global bounds: %r' % enforce_bounds)
        simulation.log.info('        Enforcing BCs: %r' % enforce_bcs)
        simulation.log.info('        Using weak BCs: %r' % use_weak_bcs)

        # Store input
        self.phi_name = phi_name
        self.phi = phi
        self.degree = degree
        self.mesh = mesh
        self.boundary_conditions = boundary_conditions
        self.use_cpp = use_cpp
        self.enforce_global_bounds = enforce_bounds
        self.enforce_boundary_conditions = enforce_bcs
        self.use_weak_bcs = use_weak_bcs
        self.num_cells_owned = mesh.topology().ghost_offset(tdim)
        self.ndim = gdim

        # No limiter needed for piecewice constant functions
        if self.degree == 0:
            self.additional_plot_funcs = []
            return

        # Alpha factors are secondary outputs
        V0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.alpha_funcs = []
        for i in range(degree):
            func = df.Function(V0)
            name = 'SlopeLimiterAlpha%d_%s' % (i + 1, phi_name)
            func.rename(name, name)
            self.alpha_funcs.append(func)
        self.additional_plot_funcs = self.alpha_funcs

        # Intermediate DG Taylor function space
        self.taylor = df.Function(V)
        self.taylor_old = df.Function(V)

        # Remove given cells from limiter
        self.limit_cell = numpy.ones(self.num_cells_owned, numpy.intc)
        if skip_cells is not None:
            for cid in skip_cells:
                self.limit_cell[cid] = 0

        self.input = SlopeLimiterInput(
            mesh, V, V0, use_cpp=use_cpp, trust_robin_dval=trust_robin_dval
        )
        if use_cpp:
            self.cpp_mod = self.input.get_cpp_mod()

    def run(self, use_weak_bcs=None):
        """
        Perform slope limiting of DG Lagrange functions
        """
        # No limiter needed for piecewice constant functions
        if self.degree == 0:
            return
        timer = df.Timer('Ocellaris HierarchalTaylorSlopeLimiter')

        # Update the Taylor function with the current Lagrange values
        lagrange_to_taylor(self.phi, self.taylor)
        taylor_arr = get_local(self.taylor)
        alpha_arrs = [alpha.vector().get_local() for alpha in self.alpha_funcs]

        # Get global bounds, see SlopeLimiterBase.set_initial_field()
        global_min, global_max = self.global_bounds

        # Update previous field values Taylor functions
        if self.phi_old is not None:
            lagrange_to_taylor(self.phi_old, self.taylor_old)
            taylor_arr_old = get_local(self.taylor_old)
        else:
            taylor_arr_old = taylor_arr

        # Get updated boundary conditions
        weak_vals = None
        use_weak_bcs = self.use_weak_bcs if use_weak_bcs is None else use_weak_bcs
        if use_weak_bcs:
            weak_vals = self.phi.vector().get_local()
        boundary_dof_type, boundary_dof_value = self.boundary_conditions.get_bcs(weak_vals)

        # Run the limiter implementation
        if self.use_cpp:
            self._run_cpp(
                taylor_arr,
                taylor_arr_old,
                alpha_arrs,
                global_min,
                global_max,
                boundary_dof_type,
                boundary_dof_value,
            )
        elif self.degree == 1 and self.ndim == 2:
            self._run_dg1(
                taylor_arr,
                taylor_arr_old,
                alpha_arrs[0],
                global_min,
                global_max,
                boundary_dof_type,
                boundary_dof_value,
            )
        elif self.degree == 2 and self.ndim == 2:
            self._run_dg2(
                taylor_arr,
                taylor_arr_old,
                alpha_arrs[0],
                alpha_arrs[1],
                global_min,
                global_max,
                boundary_dof_type,
                boundary_dof_value,
            )
        else:
            raise OcellarisError(
                'Unsupported dimension for Python version of the HierarchalTaylor limiter',
                'Only 2D is supported',
            )

        # Update the Lagrange function with the limited Taylor values
        set_local(self.taylor, taylor_arr, apply='insert')
        taylor_to_lagrange(self.taylor, self.phi)

        # Enforce boundary conditions
        if self.enforce_boundary_conditions:
            has_dbc = boundary_dof_type == self.boundary_conditions.BC_TYPE_DIRICHLET
            vals = self.phi.vector().get_local()
            vals[has_dbc] = boundary_dof_value[has_dbc]
            self.phi.vector().set_local(vals)
            self.phi.vector().apply('insert')

        # Update the secondary output arrays, alphas
        for alpha, alpha_arr in zip(self.alpha_funcs, alpha_arrs):
            alpha.vector().set_local(alpha_arr)
            alpha.vector().apply('insert')

        timer.stop()

    def _run_cpp(
        self,
        taylor_arr,
        taylor_arr_old,
        alpha_arrs,
        global_min,
        global_max,
        boundary_dof_type,
        boundary_dof_value,
    ):
        """
        Run the C++ implementation of the HierarchicalTaylor slope limiter

        The C++ versions are probably the best tested since they are fastest
        and hence most used
        """
        funcs = {
            (2, 1): self.cpp_mod.hierarchical_taylor_slope_limiter_dg1_2D,
            (2, 2): self.cpp_mod.hierarchical_taylor_slope_limiter_dg2_2D,
            (3, 1): self.cpp_mod.hierarchical_taylor_slope_limiter_dg1_3D,
            (3, 2): self.cpp_mod.hierarchical_taylor_slope_limiter_dg2_3D,
        }
        key = (self.ndim, self.degree)
        if key not in funcs:
            raise OcellarisError(
                'Unsupported dimension %d with degree %d' % key,
                'Not supported in C++ version of the HierarchalTaylor limiter',
            )

        # Update C++ input
        inp = self.input
        inp.set_global_bounds(global_min, global_max)
        inp.set_limit_cell(self.limit_cell)
        inp.set_boundary_values(
            boundary_dof_type, boundary_dof_value, self.enforce_boundary_conditions
        )

        limiter = funcs[key]
        limiter(inp.cpp_obj, taylor_arr, taylor_arr_old, *alpha_arrs)

    def _run_dg1(
        self,
        taylor_arr,
        taylor_arr_old,
        alpha_arr,
        global_min,
        global_max,
        boundary_dof_type,
        boundary_dof_value,
    ):
        """
        Perform slope limiting of a DG1 function
        """
        inp = self.input
        lagrange_arr = get_local(self.phi)
        for icell in range(self.num_cells_owned):
            dofs = inp.cell_dofs_V[icell]
            center_value = taylor_arr[dofs[0]]
            skip_this_cell = self.limit_cell[icell] == 0

            # Find the minimum slope limiter coefficient alpha
            alpha = 1.0
            if not skip_this_cell:
                for i in range(3):
                    dof = dofs[i]
                    nn = inp.num_neighbours[dof]
                    if nn == 0:
                        skip_this_cell = True
                        break

                    # Find vertex neighbours minimum and maximum values
                    minval = maxval = center_value
                    for nb in inp.neighbours[dof, :nn]:
                        nb_center_val_dof = inp.cell_dofs_V[nb][0]
                        nb_val = taylor_arr[nb_center_val_dof]
                        minval = min(minval, nb_val)
                        maxval = max(maxval, nb_val)

                        nb_val = taylor_arr_old[nb_center_val_dof]
                        minval = min(minval, nb_val)
                        maxval = max(maxval, nb_val)

                    # Modify local bounds to incorporate the global bounds
                    minval = max(minval, global_min)
                    maxval = min(maxval, global_max)
                    center_value = max(center_value, global_min)
                    center_value = min(center_value, global_max)

                    vertex_value = lagrange_arr[dof]
                    if vertex_value > center_value:
                        alpha = min(alpha, (maxval - center_value) / (vertex_value - center_value))
                    elif vertex_value < center_value:
                        alpha = min(alpha, (minval - center_value) / (vertex_value - center_value))

            if skip_this_cell:
                alpha = 1.0

            alpha_arr[inp.cell_dofs_V0[icell]] = alpha
            taylor_arr[dofs[0]] = center_value
            taylor_arr[dofs[1]] *= alpha
            taylor_arr[dofs[2]] *= alpha

    def _run_dg2(
        self,
        taylor_arr,
        taylor_arr_old,
        alpha1_arr,
        alpha2_arr,
        global_min,
        global_max,
        boundary_dof_type,
        boundary_dof_value,
    ):
        """
        Perform slope limiting of a DG2 function
        """
        # Slope limit one cell at a time
        for icell in range(self.num_cells_owned):
            dofs = self.cell_dofs_V[icell]
            assert len(dofs) == 6
            center_values = [taylor_arr[dof] for dof in dofs]
            (
                center_phi,
                center_phix,
                center_phiy,
                center_phixx,
                center_phiyy,
                center_phixy,
            ) = center_values
            skip_this_cell = self.limit_cell[icell] == 0

            cell_vertices = [self.vertex_coordinates[iv] for iv in self.vertices[icell]]
            center_pos_x = (cell_vertices[0][0] + cell_vertices[1][0] + cell_vertices[2][0]) / 3
            center_pos_y = (cell_vertices[0][1] + cell_vertices[1][1] + cell_vertices[2][1]) / 3
            assert len(cell_vertices) == 3

            # Find the minimum slope limiter coefficient alpha of the φ, dφdx and dφ/dy terms
            alpha = [1.0] * 3
            for taylor_dof in (0, 1, 2):
                if skip_this_cell:
                    break
                for ivert in range(3):
                    dof = dofs[ivert]
                    dx = cell_vertices[ivert][0] - center_pos_x
                    dy = cell_vertices[ivert][1] - center_pos_y

                    nn = self.num_neighbours[dof]
                    if nn == 0:
                        skip_this_cell = True
                        break

                    # Find vertex neighbours minimum and maximum values
                    base_value = center_values[taylor_dof]
                    minval = maxval = base_value
                    for nb in self.neighbours[dof]:
                        nb_center_val_dof = self.cell_dofs_V[nb][taylor_dof]
                        nb_val = taylor_arr[nb_center_val_dof]
                        minval = min(minval, nb_val)
                        maxval = max(maxval, nb_val)

                        nb_val = taylor_arr_old[nb_center_val_dof]
                        minval = min(minval, nb_val)
                        maxval = max(maxval, nb_val)

                    # Compute vertex value
                    if taylor_dof == 0:
                        # Modify local bounds to incorporate the global bounds
                        minval = max(minval, global_min)
                        maxval = min(maxval, global_max)
                        center_phi = max(center_phi, global_min)
                        center_phi = min(center_phi, global_max)
                        # Function value at the vertex (linear reconstruction)
                        vertex_value = center_phi + center_phix * dx + center_phiy * dy
                    elif taylor_dof == 1:
                        # Derivative in x direction at the vertex  (linear reconstruction)
                        vertex_value = center_phix + center_phixx * dx + center_phixy * dy
                    else:
                        # Derivative in y direction at the vertex  (linear reconstruction)
                        vertex_value = center_phiy + center_phiyy * dy + center_phixy * dx

                    # Compute the slope limiter coefficient alpha
                    if vertex_value > base_value:
                        a = (maxval - base_value) / (vertex_value - base_value)
                    elif vertex_value < base_value:
                        a = (minval - base_value) / (vertex_value - base_value)
                    else:
                        a = 1
                    alpha[taylor_dof] = min(alpha[taylor_dof], a)

            if skip_this_cell:
                alpha1 = alpha2 = 1.0
            else:
                alpha2 = min(alpha[1], alpha[2])
                alpha1 = max(alpha[0], alpha2)

            taylor_arr[dofs[0]] = center_phi
            taylor_arr[dofs[1]] *= alpha1
            taylor_arr[dofs[2]] *= alpha1
            taylor_arr[dofs[3]] *= alpha2
            taylor_arr[dofs[4]] *= alpha2
            taylor_arr[dofs[5]] *= alpha2

            dof_dg0 = self.cell_dofs_V0[icell]
            alpha1_arr[dof_dg0] = alpha1
            alpha2_arr[dof_dg0] = alpha2
