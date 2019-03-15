# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from ocellaris.utils import lagrange_to_taylor, taylor_to_lagrange
import pytest


@pytest.mark.parametrize("degree", [1, 2])
def test_taylor_projections_2D(degree):
    """
    Check that the 2D Lagrange -> Taylor -> Lagrange projection is an identity projection
    """
    mesh = dolfin.UnitSquareMesh(4, 4)
    V = dolfin.FunctionSpace(mesh, 'DG', degree)

    # Setup Lagrange function with random data
    f1 = dolfin.Function(V)
    N = f1.vector().local_size()
    f1.vector().set_local(numpy.random.rand(N))

    # Convert to Taylor representation
    f2 = dolfin.Function(V)
    lagrange_to_taylor(f1, f2)

    # Create a new Lagrange function from the Taylor representation
    f3 = dolfin.Function(V)
    taylor_to_lagrange(f2, f3)

    error = dolfin.errornorm(f1, f3, degree_rise=0)
    assert error < 1e-15


@pytest.mark.parametrize("degree", [1, 2])
def test_taylor_projections_3D(degree):
    """
    Check that the 3D Lagrange -> Taylor -> Lagrange projection is an identity projection
    """
    mesh = dolfin.UnitCubeMesh(2, 2, 2)
    V = dolfin.FunctionSpace(mesh, 'DG', degree)

    # Setup Lagrange function with random data
    f1 = dolfin.Function(V)
    N = f1.vector().local_size()
    f1.vector().set_local(numpy.random.rand(N))

    # Convert to Taylor representation
    f2 = dolfin.Function(V)
    lagrange_to_taylor(f1, f2)

    # Create a new Lagrange function from the Taylor representation
    f3 = dolfin.Function(V)
    taylor_to_lagrange(f2, f3)

    error = dolfin.errornorm(f1, f3, degree_rise=0)
    assert error < 1e-15


def make_taylor_func(u, coeffs):
    """
    Modify the Lagrange function u such that for every cell
    the function is
        
        A, B, C, ... = coeffs
        f_cell = (Â + B*(x-xc) + C*(y-yc) + D*(x-xc)**2/2 + 
                  E*(y-yc)**2/2 + F*(x-xc)*(y-yc))
    
    Where Â is such that the cell average is A 
    """
    V = u.function_space()
    dm = V.dofmap()
    mesh = V.mesh()
    gdim = mesh.geometry().dim()
    dofs_x = V.tabulate_dof_coordinates().reshape((-1, gdim))
    vals = u.vector().get_local()

    for cell in dolfin.cells(mesh):
        cell_dofs = dm.cell_dofs(cell.index())
        N = len(cell_dofs)

        # Find the cell center
        if N in (3, 6):
            # Triangle cell center
            d0, d1, d2 = cell_dofs[:3]
            xc = (dofs_x[d0] + dofs_x[d1] + dofs_x[d2]) / 3
        elif N in (4, 10):
            # Tetrahedron cell center
            d0, d1, d2, d3 = cell_dofs[:4]
            xc = (dofs_x[d0] + dofs_x[d1] + dofs_x[d2] + dofs_x[d3]) / 4

        # Assign the non-const part
        for dof in cell_dofs:
            pos = dofs_x[dof]
            if N == 3:
                # Linear triangle
                A, B, C = coeffs
                vals[dof] = B * (pos[0] - xc[0]) + C * (pos[1] - xc[1])
            elif N == 6:
                # Quadratic triangle
                A, B, C, D, E, F = coeffs
                vals[dof] = B * (pos[0] - xc[0]) + C * (pos[1] - xc[1])
                vals[dof] += D * (pos[0] - xc[0]) ** 2 / 2
                vals[dof] += E * (pos[1] - xc[1]) ** 2 / 2
                vals[dof] += F * (pos[0] - xc[0]) * (pos[1] - xc[1])
            elif N == 4:
                # Linear tetrahedron
                A, B, C, D = coeffs
                vals[dof] = (
                    B * (pos[0] - xc[0]) + C * (pos[1] - xc[1]) + D * (pos[2] - xc[2])
                )
            elif N == 10:
                # Quadratic tetrahedron
                A, B, C, D, E, F, G, H, I, J = coeffs
                vals[dof] = (
                    B * (pos[0] - xc[0]) + C * (pos[1] - xc[1]) + D * (pos[2] - xc[2])
                )
                vals[dof] += E * (pos[0] - xc[0]) ** 2 / 2
                vals[dof] += F * (pos[1] - xc[1]) ** 2 / 2
                vals[dof] += G * (pos[2] - xc[2]) ** 2 / 2
                vals[dof] += H * (pos[0] - xc[0]) * (pos[1] - xc[1])
                vals[dof] += I * (pos[0] - xc[0]) * (pos[2] - xc[2])
                vals[dof] += J * (pos[1] - xc[1]) * (pos[2] - xc[2])

        # Assign the const part
        if N == 3:
            avg_weights = [1 / 3, 1 / 3, 1 / 3]
        elif N == 6:
            avg_weights = [0, 0, 0, 1 / 3, 1 / 3, 1 / 3]
        elif N == 4:
            avg_weights = [1 / 4, 1 / 4, 1 / 4, 1 / 4]
        elif N == 10:
            avg_weights = [
                -1 / 20,
                -1 / 20,
                -1 / 20,
                -1 / 20,
                1 / 5,
                1 / 5,
                1 / 5,
                1 / 5,
                1 / 5,
                1 / 5,
            ]
        assert abs(1 - sum(avg_weights)) < 1e-14
        assert len(avg_weights) == N
        avg = numpy.dot(avg_weights, vals[cell_dofs])
        for dof in cell_dofs:
            vals[dof] += A - avg

    u.vector().set_local(vals)
    u.vector().apply('insert')


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("degree", [1, 2])
def test_taylor_values(dim, degree):
    """
    Check that the Lagrange -> Taylor projection gives the correct Taylor values
    """
    if dim == 2:
        mesh = dolfin.UnitSquareMesh(4, 4)
    else:
        mesh = dolfin.UnitCubeMesh(2, 2, 2)

    # Setup Lagrange function with given derivatives and constants
    V = dolfin.FunctionSpace(mesh, 'DG', degree)
    u = dolfin.Function(V)
    if dim == 2:
        coeffs = [1, 2, -3.0] if degree == 1 else [1, 2, -3.0, -2.5, 4.2, -1.0]
    else:
        coeffs = (
            [1, 2, -3.0, 2.5]
            if degree == 1
            else [1, 2, -3.0, 2.5, -1.3, 4.2, -1.0, -4.2, 2.66, 3.14]
        )
    make_taylor_func(u, coeffs)

    # Convert to Taylor
    t = dolfin.Function(V)
    lagrange_to_taylor(u, t)

    # Check that the target values are obtained
    dm = V.dofmap()
    vals = t.vector().get_local()
    for cell in dolfin.cells(mesh):
        cell_dofs = dm.cell_dofs(cell.index())
        cell_vals = vals[cell_dofs]
        assert all(abs(cell_vals - coeffs) < 1e-13)
