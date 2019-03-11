# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
Convert to and from a Taylor basis

This code follows the definition of cell average Taylor DG polynomial expansion
as described in Dmitri Kuzmin (2010) "A vertex-based hierarchial slope limiter
for p-adaptive discontinuous Galerkin methods"

Code generation the transformation is implemented using SymPy and barycentric
coordinates and can be found in documentation/notebooks/barycentric.ipynb
"""
import numpy
from dolfin import cells
from ocellaris.utils import ocellaris_error


CACHE = {}


def lagrange_to_taylor(u, t):
    """
    Convert a Lagrange function space function into a DG Taylor
    function space function. The results are stored in t which
    should be a function space with an appropriate number of dofs
    per cell (t in DG2 if u in DG2 etc)

    Note: Most possible dolfin operations on the function t will
    be wrong since dolfin does not support Taylor basis DG elements.
    The t function should hence be read and modified by specially
    crafted code only!
    """
    V = u.function_space()
    mesh = V.mesh()
    degree = V.ufl_element().degree()
    ndim = mesh.geometry().dim()

    if 'mesh_hash' not in CACHE or CACHE['mesh_hash'] != mesh.hash():
        CACHE.clear()
        CACHE['mesh'] = mesh
        CACHE['mesh_hash'] = mesh.hash()

    # Get cached conversion matrices
    key = ('lagrange_to_taylor_matrices', degree)
    if key not in CACHE:
        if degree == 1 and ndim == 2:
            CACHE[key] = DG1_to_taylor_matrix_2D(V)
        elif degree == 1 and ndim == 3:
            CACHE[key] = DG1_to_taylor_matrix_3D(V)
        elif degree == 2 and ndim == 2:
            CACHE[key] = DG2_to_taylor_matrix_2D(V)
        elif degree == 2 and ndim == 3:
            CACHE[key] = DG2_to_taylor_matrix_3D(V)
        else:
            ocellaris_error(
                'DG Lagrange to DG Taylor converter error',
                'Polynomial degree %d not supported' % degree,
            )
    lagrange_to_taylor_matrices = CACHE[key]
    Ncells = lagrange_to_taylor_matrices.shape[0]

    # Get cached cell dofs
    key2 = ('cell_dofs', degree)
    if key2 not in CACHE:
        dm = V.dofmap()
        cell_dofs = [dm.cell_dofs(i) for i in range(Ncells)]
        CACHE[key2] = numpy.array(cell_dofs, int)
    cell_dofs = CACHE[key2]

    # Apply the conversion matrices for all cells by use of the stacked dot
    # behaviour of matmul (slightly faster than einsum 'ijk,ik->ij')
    all_vals_lagrange = u.vector().get_local()
    lagrange_vectors = all_vals_lagrange.take(cell_dofs)
    res = numpy.matmul(
        lagrange_to_taylor_matrices, lagrange_vectors[:, :, None]
    ).squeeze()

    # Put the results into the right indices in the Taylor function's vector
    all_vals_taylor = numpy.zeros_like(all_vals_lagrange)
    all_vals_taylor[cell_dofs] = res
    t.vector().set_local(all_vals_taylor)
    t.vector().apply('insert')


def taylor_to_lagrange(t, u):
    """
    Convert a DG Taylor function space function into a Lagrange
    function space function. The results are stored in u which
    should be a function space with an appropriate number of dofs
    per cell (u in DG2 if t in DG2 etc)
    """
    V = u.function_space()
    mesh = V.mesh()
    degree = V.ufl_element().degree()
    ndim = mesh.geometry().dim()

    if 'mesh_hash' not in CACHE or CACHE['mesh_hash'] != mesh.hash():
        CACHE.clear()
        CACHE['mesh'] = mesh
        CACHE['mesh_hash'] = mesh.hash()

    # Get cached conversion matrices
    key = ('taylor_to_lagrange_matrices', degree)
    if key not in CACHE:
        if degree == 1 and ndim == 2:
            CACHE[key] = taylor_to_DG1_matrix_2D(V)
        elif degree == 1 and ndim == 3:
            CACHE[key] = taylor_to_DG1_matrix_3D(V)
        elif degree == 2 and ndim == 2:
            CACHE[key] = taylor_to_DG2_matrix_2D(V)
        elif degree == 2 and ndim == 3:
            CACHE[key] = taylor_to_DG2_matrix_3D(V)
        else:
            ocellaris_error(
                'DG Taylor to DG Lagrange converter error',
                'Polynomial degree %d not supported' % degree,
            )
    taylor_to_lagrange_matrices = CACHE[key]
    Ncells = taylor_to_lagrange_matrices.shape[0]

    # Get cached cell dofs
    key2 = ('cell_dofs', degree)
    if key2 not in CACHE:
        dm = V.dofmap()
        cell_dofs = [dm.cell_dofs(i) for i in range(Ncells)]
        CACHE[key2] = numpy.array(cell_dofs, int)
    cell_dofs = CACHE[key2]

    # Apply the conversion matrices for all cells by use of the stacked dot
    # behaviour of matmul (slightly faster than einsum 'ijk,ik->ij')
    all_vals_taylor = t.vector().get_local()
    taylor_vectors = all_vals_taylor.take(cell_dofs)
    res = numpy.matmul(
        taylor_to_lagrange_matrices, taylor_vectors[:, :, None]
    ).squeeze()

    # Put the results into the right indices in the Taylor function's vector
    all_vals_lagrange = numpy.zeros_like(all_vals_taylor)
    all_vals_lagrange[cell_dofs] = res
    u.vector().set_local(all_vals_lagrange)
    u.vector().apply('insert')


##########################################################################
# DG Lagrange to Taylor


def DG1_to_taylor_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG1 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    A = numpy.zeros((num_cells_owned, 3, 3), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x = numpy.zeros((3, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]

        ###############################
        # From sympy code gen, see documentation/notebooks/barycentric.ipynb

        ((x1, y1), (x2, y2), (x3, y3)) = x
        D = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)

        # Value at xc (also the cell average value)
        A[icell, 0, 0] = 1 / 3
        A[icell, 0, 1] = 1 / 3
        A[icell, 0, 2] = 1 / 3

        # d/dx
        A[icell, 1, 0] = (y2 - y3) / D
        A[icell, 1, 1] = (-y1 + y3) / D
        A[icell, 1, 2] = (y1 - y2) / D

        # d/dy
        A[icell, 2, 0] = (-x2 + x3) / D
        A[icell, 2, 1] = (x1 - x3) / D
        A[icell, 2, 2] = (-x1 + x2) / D

    return A


def DG1_to_taylor_matrix_3D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG1 function space V in 3D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    A = numpy.zeros((num_cells_owned, 4, 4), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x = numpy.zeros((4, 3), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = vertices[verts[3]]

        ###############################
        # From sympy code gen, see documentation/notebooks/barycentric.ipynb

        ((x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4)) = x
        F = (
            x1 * y2 * z3
            - x1 * y2 * z4
            - x1 * y3 * z2
            + x1 * y3 * z4
            + x1 * y4 * z2
            - x1 * y4 * z3
            - x2 * y1 * z3
            + x2 * y1 * z4
            + x2 * y3 * z1
            - x2 * y3 * z4
            - x2 * y4 * z1
            + x2 * y4 * z3
            + x3 * y1 * z2
            - x3 * y1 * z4
            - x3 * y2 * z1
            + x3 * y2 * z4
            + x3 * y4 * z1
            - x3 * y4 * z2
            - x4 * y1 * z2
            + x4 * y1 * z3
            + x4 * y2 * z1
            - x4 * y2 * z3
            - x4 * y3 * z1
            + x4 * y3 * z2
        )

        # Value at xc (also the cell average value)
        A[icell, 0, 0] = 1 / 4
        A[icell, 0, 1] = 1 / 4
        A[icell, 0, 2] = 1 / 4
        A[icell, 0, 3] = 1 / 4

        # d/dx
        A[icell, 1, 0] = (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3) / F
        A[icell, 1, 1] = (
            -y1 * z3 + y1 * z4 + y3 * z1 - y3 * z4 - y4 * z1 + y4 * z3
        ) / F
        A[icell, 1, 2] = (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2) / F
        A[icell, 1, 3] = (
            -y1 * z2 + y1 * z3 + y2 * z1 - y2 * z3 - y3 * z1 + y3 * z2
        ) / F

        # d/dy
        A[icell, 2, 0] = (
            -x2 * z3 + x2 * z4 + x3 * z2 - x3 * z4 - x4 * z2 + x4 * z3
        ) / F
        A[icell, 2, 1] = (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3) / F
        A[icell, 2, 2] = (
            -x1 * z2 + x1 * z4 + x2 * z1 - x2 * z4 - x4 * z1 + x4 * z2
        ) / F
        A[icell, 2, 3] = (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) / F

        # d/dz
        A[icell, 3, 0] = (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3) / F
        A[icell, 3, 1] = (
            -x1 * y3 + x1 * y4 + x3 * y1 - x3 * y4 - x4 * y1 + x4 * y3
        ) / F
        A[icell, 3, 2] = (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2) / F
        A[icell, 3, 3] = (
            -x1 * y2 + x1 * y3 + x2 * y1 - x2 * y3 - x3 * y1 + x3 * y2
        ) / F

    return A


def DG2_to_taylor_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG2 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    A = numpy.zeros((num_cells_owned, 6, 6), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x = numpy.zeros((3, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]

        ###############################
        # From sympy code gen, see documentation/notebooks/barycentric.ipynb

        ((x1, y1), (x2, y2), (x3, y3)) = x
        D = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)

        # Cell average value
        A[icell, 0, 3] = 1 / 3
        A[icell, 0, 4] = 1 / 3
        A[icell, 0, 5] = 1 / 3

        # d/dx
        A[icell, 1, 0] = (y2 - y3) / (3 * D)
        A[icell, 1, 1] = (-y1 + y3) / (3 * D)
        A[icell, 1, 2] = (y1 - y2) / (3 * D)
        A[icell, 1, 3] = 4 * (-y2 + y3) / (3 * D)
        A[icell, 1, 4] = 4 * (y1 - y3) / (3 * D)
        A[icell, 1, 5] = 4 * (-y1 + y2) / (3 * D)

        # d/dy
        A[icell, 2, 0] = (-x2 + x3) / (3 * D)
        A[icell, 2, 1] = (x1 - x3) / (3 * D)
        A[icell, 2, 2] = (-x1 + x2) / (3 * D)
        A[icell, 2, 3] = 4 * (x2 - x3) / (3 * D)
        A[icell, 2, 4] = 4 * (-x1 + x3) / (3 * D)
        A[icell, 2, 5] = 4 * (x1 - x2) / (3 * D)

        # d/dx^2
        A[icell, 3, 0] = 4 * (y2 - y3) ** 2 / D ** 2
        A[icell, 3, 1] = 4 * (y1 - y3) ** 2 / D ** 2
        A[icell, 3, 2] = 4 * (y1 - y2) ** 2 / D ** 2
        A[icell, 3, 3] = -8 * (y1 - y2) * (y1 - y3) / D ** 2
        A[icell, 3, 4] = 8 * (y1 - y2) * (y2 - y3) / D ** 2
        A[icell, 3, 5] = -8 * (y1 - y3) * (y2 - y3) / D ** 2

        # d/dy^2
        A[icell, 4, 0] = 4 * (x2 - x3) ** 2 / D ** 2
        A[icell, 4, 1] = 4 * (x1 - x3) ** 2 / D ** 2
        A[icell, 4, 2] = 4 * (x1 - x2) ** 2 / D ** 2
        A[icell, 4, 3] = -8 * (x1 - x2) * (x1 - x3) / D ** 2
        A[icell, 4, 4] = 8 * (x1 - x2) * (x2 - x3) / D ** 2
        A[icell, 4, 5] = -8 * (x1 - x3) * (x2 - x3) / D ** 2

        # d/dx*dy
        A[icell, 5, 0] = -4 * (x2 - x3) * (y2 - y3) / D ** 2
        A[icell, 5, 1] = -4 * (x1 - x3) * (y1 - y3) / D ** 2
        A[icell, 5, 2] = -4 * (x1 - x2) * (y1 - y2) / D ** 2
        A[icell, 5, 3] = 4 * ((x1 - x2) * (y1 - y3) + (x1 - x3) * (y1 - y2)) / D ** 2
        A[icell, 5, 4] = (
            -(4 * (x1 - x2) * (y2 - y3) + 4 * (x2 - x3) * (y1 - y2)) / D ** 2
        )
        A[icell, 5, 5] = 4 * ((x1 - x3) * (y2 - y3) + (x2 - x3) * (y1 - y3)) / D ** 2

    return A


def DG2_to_taylor_matrix_3D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG1 function space V in 3D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    A = numpy.zeros((num_cells_owned, 10, 10), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x = numpy.zeros((4, 3), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = vertices[verts[3]]

        ###############################
        # From sympy code gen, see documentation/notebooks/barycentric.ipynb

        ((x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4)) = x
        F = (
            x1 * y2 * z3
            - x1 * y2 * z4
            - x1 * y3 * z2
            + x1 * y3 * z4
            + x1 * y4 * z2
            - x1 * y4 * z3
            - x2 * y1 * z3
            + x2 * y1 * z4
            + x2 * y3 * z1
            - x2 * y3 * z4
            - x2 * y4 * z1
            + x2 * y4 * z3
            + x3 * y1 * z2
            - x3 * y1 * z4
            - x3 * y2 * z1
            + x3 * y2 * z4
            + x3 * y4 * z1
            - x3 * y4 * z2
            - x4 * y1 * z2
            + x4 * y1 * z3
            + x4 * y2 * z1
            - x4 * y2 * z3
            - x4 * y3 * z1
            + x4 * y3 * z2
        )

        # Cell average value
        A[icell, 0, 0] = -1 / 20
        A[icell, 0, 1] = -1 / 20
        A[icell, 0, 2] = -1 / 20
        A[icell, 0, 3] = -1 / 20
        A[icell, 0, 4] = 1 / 5
        A[icell, 0, 5] = 1 / 5
        A[icell, 0, 6] = 1 / 5
        A[icell, 0, 7] = 1 / 5
        A[icell, 0, 8] = 1 / 5
        A[icell, 0, 9] = 1 / 5

        # d/dx
        A[icell, 1, 0] = 0
        A[icell, 1, 1] = 0
        A[icell, 1, 2] = 0
        A[icell, 1, 3] = 0
        A[icell, 1, 4] = (
            y1 * z3
            - y1 * z4
            - y2 * z3
            + y2 * z4
            - y3 * z1
            + y3 * z2
            + y4 * z1
            - y4 * z2
        ) / F
        A[icell, 1, 5] = (
            -y1 * z2
            + y1 * z4
            + y2 * z1
            - y2 * z3
            + y3 * z2
            - y3 * z4
            - y4 * z1
            + y4 * z3
        ) / F
        A[icell, 1, 6] = (
            y1 * z2
            - y1 * z3
            - y2 * z1
            + y2 * z4
            + y3 * z1
            - y3 * z4
            - y4 * z2
            + y4 * z3
        ) / F
        A[icell, 1, 7] = (
            -y1 * z2
            + y1 * z3
            + y2 * z1
            - y2 * z4
            - y3 * z1
            + y3 * z4
            + y4 * z2
            - y4 * z3
        ) / F
        A[icell, 1, 8] = (
            y1 * z2
            - y1 * z4
            - y2 * z1
            + y2 * z3
            - y3 * z2
            + y3 * z4
            + y4 * z1
            - y4 * z3
        ) / F
        A[icell, 1, 9] = (
            -y1 * z3
            + y1 * z4
            + y2 * z3
            - y2 * z4
            + y3 * z1
            - y3 * z2
            - y4 * z1
            + y4 * z2
        ) / F

        # d/dy
        A[icell, 2, 0] = 0
        A[icell, 2, 1] = 0
        A[icell, 2, 2] = 0
        A[icell, 2, 3] = 0
        A[icell, 2, 4] = (
            -x1 * z3
            + x1 * z4
            + x2 * z3
            - x2 * z4
            + x3 * z1
            - x3 * z2
            - x4 * z1
            + x4 * z2
        ) / F
        A[icell, 2, 5] = (
            x1 * z2
            - x1 * z4
            - x2 * z1
            + x2 * z3
            - x3 * z2
            + x3 * z4
            + x4 * z1
            - x4 * z3
        ) / F
        A[icell, 2, 6] = (
            -x1 * z2
            + x1 * z3
            + x2 * z1
            - x2 * z4
            - x3 * z1
            + x3 * z4
            + x4 * z2
            - x4 * z3
        ) / F
        A[icell, 2, 7] = (
            x1 * z2
            - x1 * z3
            - x2 * z1
            + x2 * z4
            + x3 * z1
            - x3 * z4
            - x4 * z2
            + x4 * z3
        ) / F
        A[icell, 2, 8] = (
            -x1 * z2
            + x1 * z4
            + x2 * z1
            - x2 * z3
            + x3 * z2
            - x3 * z4
            - x4 * z1
            + x4 * z3
        ) / F
        A[icell, 2, 9] = (
            x1 * z3
            - x1 * z4
            - x2 * z3
            + x2 * z4
            - x3 * z1
            + x3 * z2
            + x4 * z1
            - x4 * z2
        ) / F

        # d/dz
        A[icell, 3, 0] = 0
        A[icell, 3, 1] = 0
        A[icell, 3, 2] = 0
        A[icell, 3, 3] = 0
        A[icell, 3, 4] = (
            x1 * y3
            - x1 * y4
            - x2 * y3
            + x2 * y4
            - x3 * y1
            + x3 * y2
            + x4 * y1
            - x4 * y2
        ) / F
        A[icell, 3, 5] = (
            -x1 * y2
            + x1 * y4
            + x2 * y1
            - x2 * y3
            + x3 * y2
            - x3 * y4
            - x4 * y1
            + x4 * y3
        ) / F
        A[icell, 3, 6] = (
            x1 * y2
            - x1 * y3
            - x2 * y1
            + x2 * y4
            + x3 * y1
            - x3 * y4
            - x4 * y2
            + x4 * y3
        ) / F
        A[icell, 3, 7] = (
            -x1 * y2
            + x1 * y3
            + x2 * y1
            - x2 * y4
            - x3 * y1
            + x3 * y4
            + x4 * y2
            - x4 * y3
        ) / F
        A[icell, 3, 8] = (
            x1 * y2
            - x1 * y4
            - x2 * y1
            + x2 * y3
            - x3 * y2
            + x3 * y4
            + x4 * y1
            - x4 * y3
        ) / F
        A[icell, 3, 9] = (
            -x1 * y3
            + x1 * y4
            + x2 * y3
            - x2 * y4
            + x3 * y1
            - x3 * y2
            - x4 * y1
            + x4 * y2
        ) / F

        # d/dx^2
        A[icell, 4, 0] = (
            4
            * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3) ** 2
            / F ** 2
        )
        A[icell, 4, 1] = (
            4
            * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3) ** 2
            / F ** 2
        )
        A[icell, 4, 2] = (
            4
            * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2) ** 2
            / F ** 2
        )
        A[icell, 4, 3] = (
            4
            * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2) ** 2
            / F ** 2
        )
        A[icell, 4, 4] = (
            -8
            * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            / F ** 2
        )
        A[icell, 4, 5] = (
            8
            * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            / F ** 2
        )
        A[icell, 4, 6] = (
            -8
            * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            / F ** 2
        )
        A[icell, 4, 7] = (
            -8
            * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
            / F ** 2
        )
        A[icell, 4, 8] = (
            8
            * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
            / F ** 2
        )
        A[icell, 4, 9] = (
            -8
            * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
            / F ** 2
        )

        # d/dy^2
        A[icell, 5, 0] = (
            4
            * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3) ** 2
            / F ** 2
        )
        A[icell, 5, 1] = (
            4
            * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3) ** 2
            / F ** 2
        )
        A[icell, 5, 2] = (
            4
            * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2) ** 2
            / F ** 2
        )
        A[icell, 5, 3] = (
            4
            * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2) ** 2
            / F ** 2
        )
        A[icell, 5, 4] = (
            -8
            * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
            / F ** 2
        )
        A[icell, 5, 5] = (
            8
            * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
            / F ** 2
        )
        A[icell, 5, 6] = (
            -8
            * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
            * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
            / F ** 2
        )
        A[icell, 5, 7] = (
            -8
            * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
            / F ** 2
        )
        A[icell, 5, 8] = (
            8
            * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
            * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
            / F ** 2
        )
        A[icell, 5, 9] = (
            -8
            * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
            * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
            / F ** 2
        )

        # d/dz^2
        A[icell, 6, 0] = (
            4
            * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3) ** 2
            / F ** 2
        )
        A[icell, 6, 1] = (
            4
            * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3) ** 2
            / F ** 2
        )
        A[icell, 6, 2] = (
            4
            * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2) ** 2
            / F ** 2
        )
        A[icell, 6, 3] = (
            4
            * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) ** 2
            / F ** 2
        )
        A[icell, 6, 4] = (
            -8
            * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
            * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
            / F ** 2
        )
        A[icell, 6, 5] = (
            8
            * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
            * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
            / F ** 2
        )
        A[icell, 6, 6] = (
            -8
            * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
            * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
            / F ** 2
        )
        A[icell, 6, 7] = (
            -8
            * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
            * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            / F ** 2
        )
        A[icell, 6, 8] = (
            8
            * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
            * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            / F ** 2
        )
        A[icell, 6, 9] = (
            -8
            * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
            * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            / F ** 2
        )

        # d/dx*dy
        A[icell, 7, 0] = (
            -4
            * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
            * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
            / F ** 2
        )
        A[icell, 7, 1] = (
            -4
            * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
            * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            / F ** 2
        )
        A[icell, 7, 2] = (
            -4
            * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
            * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            / F ** 2
        )
        A[icell, 7, 3] = (
            -4
            * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            / F ** 2
        )
        A[icell, 7, 4] = (
            4
            * (
                (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
                * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
                + (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
                * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            )
            / F ** 2
        )
        A[icell, 7, 5] = (
            -(
                4
                * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
                * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
                + 4
                * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
                * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            )
            / F ** 2
        )
        A[icell, 7, 6] = (
            4
            * (
                (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
                * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
                + (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
                * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            )
            / F ** 2
        )
        A[icell, 7, 7] = (
            4
            * (
                (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
                * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
                + (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
                * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            )
            / F ** 2
        )
        A[icell, 7, 8] = (
            -(
                4
                * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
                * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
                + 4
                * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
                * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            )
            / F ** 2
        )
        A[icell, 7, 9] = (
            4
            * (
                (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
                * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
                + (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
                * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            )
            / F ** 2
        )

        # d/dx*dz
        A[icell, 8, 0] = (
            4
            * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
            / F ** 2
        )
        A[icell, 8, 1] = (
            4
            * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
            * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            / F ** 2
        )
        A[icell, 8, 2] = (
            4
            * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
            * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            / F ** 2
        )
        A[icell, 8, 3] = (
            4
            * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
            * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            / F ** 2
        )
        A[icell, 8, 4] = (
            -(
                4
                * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
                * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
                + 4
                * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
                * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            )
            / F ** 2
        )
        A[icell, 8, 5] = (
            4
            * (
                (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
                * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
                + (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
                * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            )
            / F ** 2
        )
        A[icell, 8, 6] = (
            -(
                4
                * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
                * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
                + 4
                * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
                * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            )
            / F ** 2
        )
        A[icell, 8, 7] = (
            -(
                4
                * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
                * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
                + 4
                * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
                * (y1 * z2 - y1 * z3 - y2 * z1 + y2 * z3 + y3 * z1 - y3 * z2)
            )
            / F ** 2
        )
        A[icell, 8, 8] = (
            4
            * (
                (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
                * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
                + (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
                * (y1 * z2 - y1 * z4 - y2 * z1 + y2 * z4 + y4 * z1 - y4 * z2)
            )
            / F ** 2
        )
        A[icell, 8, 9] = (
            -(
                4
                * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
                * (y2 * z3 - y2 * z4 - y3 * z2 + y3 * z4 + y4 * z2 - y4 * z3)
                + 4
                * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
                * (y1 * z3 - y1 * z4 - y3 * z1 + y3 * z4 + y4 * z1 - y4 * z3)
            )
            / F ** 2
        )

        # d/dy*dz
        A[icell, 9, 0] = (
            -4
            * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
            / F ** 2
        )
        A[icell, 9, 1] = (
            -4
            * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
            * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
            / F ** 2
        )
        A[icell, 9, 2] = (
            -4
            * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
            * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
            / F ** 2
        )
        A[icell, 9, 3] = (
            -4
            * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
            * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            / F ** 2
        )
        A[icell, 9, 4] = (
            4
            * (
                (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
                * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
                + (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
                * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            )
            / F ** 2
        )
        A[icell, 9, 5] = (
            -(
                4
                * (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
                * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
                + 4
                * (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
                * (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
            )
            / F ** 2
        )
        A[icell, 9, 6] = (
            4
            * (
                (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
                * (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
                + (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
                * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
            )
            / F ** 2
        )
        A[icell, 9, 7] = (
            4
            * (
                (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2)
                * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
                + (x1 * z2 - x1 * z3 - x2 * z1 + x2 * z3 + x3 * z1 - x3 * z2)
                * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            )
            / F ** 2
        )
        A[icell, 9, 8] = (
            -(
                4
                * (x1 * y2 - x1 * y4 - x2 * y1 + x2 * y4 + x4 * y1 - x4 * y2)
                * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
                + 4
                * (x1 * z2 - x1 * z4 - x2 * z1 + x2 * z4 + x4 * z1 - x4 * z2)
                * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            )
            / F ** 2
        )
        A[icell, 9, 9] = (
            4
            * (
                (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3)
                * (x2 * z3 - x2 * z4 - x3 * z2 + x3 * z4 + x4 * z2 - x4 * z3)
                + (x1 * z3 - x1 * z4 - x3 * z1 + x3 * z4 + x4 * z1 - x4 * z3)
                * (x2 * y3 - x2 * y4 - x3 * y2 + x3 * y4 + x4 * y2 - x4 * y3)
            )
            / F ** 2
        )

    return A


##########################################################################
# Taylor to DG Lagrange


def taylor_to_DG1_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG1 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    x = numpy.zeros((3, 2), float)
    A = numpy.zeros((num_cells_owned, 3, 3), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        xc = (x[0] + x[1] + x[2]) / 3

        for i in range(3):
            dx, dy = x[i, 0] - xc[0], x[i, 1] - xc[1]
            A[icell, i, 0] = 1
            A[icell, i, 1] = dx
            A[icell, i, 2] = dy

    return A


def taylor_to_DG1_matrix_3D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG1 function space V in 3D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    x = numpy.zeros((4, 3), float)
    A = numpy.zeros((num_cells_owned, 4, 4), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = vertices[verts[3]]
        xc = (x[0] + x[1] + x[2] + x[3]) / 4

        for i in range(4):
            dx, dy, dz = x[i, 0] - xc[0], x[i, 1] - xc[1], x[i, 2] - xc[2]
            A[icell, i, 0] = 1
            A[icell, i, 1] = dx
            A[icell, i, 2] = dy
            A[icell, i, 3] = dz

    return A


def taylor_to_DG2_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG2 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    x = numpy.zeros((6, 2), float)
    A = numpy.zeros((num_cells_owned, 6, 6), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = (x[1] + x[2]) / 2
        x[4] = (x[0] + x[2]) / 2
        x[5] = (x[0] + x[1]) / 2
        xc = (x[0] + x[1] + x[2]) / 3

        # Code generated by the sympy code included below
        ((x1, y1), (x2, y2), (x3, y3)) = x[:3]
        bar_xx = (
            x1 ** 2 / 36
            - x1 * x2 / 36
            - x1 * x3 / 36
            + x2 ** 2 / 36
            - x2 * x3 / 36
            + x3 ** 2 / 36
        )
        bar_yy = (
            y1 ** 2 / 36
            - y1 * y2 / 36
            - y1 * y3 / 36
            + y2 ** 2 / 36
            - y2 * y3 / 36
            + y3 ** 2 / 36
        )
        bar_xy = (
            x1 * y1 / 18
            - x1 * y2 / 36
            - x1 * y3 / 36
            - x2 * y1 / 36
            + x2 * y2 / 18
            - x2 * y3 / 36
            - x3 * y1 / 36
            - x3 * y2 / 36
            + x3 * y3 / 18
        )

        for i in range(6):
            dx, dy = x[i, 0] - xc[0], x[i, 1] - xc[1]
            A[icell, i, 0] = 1
            A[icell, i, 1] = dx
            A[icell, i, 2] = dy
            A[icell, i, 3] = dx ** 2 / 2 - bar_xx
            A[icell, i, 4] = dy ** 2 / 2 - bar_yy
            A[icell, i, 5] = dx * dy - bar_xy

    return A


def taylor_to_DG2_matrix_3D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG2 function space V in 3D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()

    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)

    x = numpy.zeros((10, 3), float)
    A = numpy.zeros((num_cells_owned, 10, 10), float)
    for cell in cells(mesh):
        icell = cell.index()

        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = vertices[verts[3]]
        x[4] = (x[2] + x[3]) / 2
        x[5] = (x[1] + x[3]) / 2
        x[6] = (x[1] + x[2]) / 2
        x[7] = (x[0] + x[3]) / 2
        x[8] = (x[0] + x[2]) / 2
        x[9] = (x[0] + x[1]) / 2
        xc = (x[0] + x[1] + x[2] + x[3]) / 4

        # Code generated by the sympy code included below
        ((x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4)) = x[:4]
        bar_xx = (
            (-3 * x1 + x2 + x3 + x4) ** 2 / 640
            + (x1 - 3 * x2 + x3 + x4) ** 2 / 640
            + (x1 + x2 - 3 * x3 + x4) ** 2 / 640
            + (x1 + x2 + x3 - 3 * x4) ** 2 / 640
        )
        bar_yy = (
            (-3 * y1 + y2 + y3 + y4) ** 2 / 640
            + (y1 - 3 * y2 + y3 + y4) ** 2 / 640
            + (y1 + y2 - 3 * y3 + y4) ** 2 / 640
            + (y1 + y2 + y3 - 3 * y4) ** 2 / 640
        )
        bar_zz = (
            (-3 * z1 + z2 + z3 + z4) ** 2 / 640
            + (z1 - 3 * z2 + z3 + z4) ** 2 / 640
            + (z1 + z2 - 3 * z3 + z4) ** 2 / 640
            + (z1 + z2 + z3 - 3 * z4) ** 2 / 640
        )
        bar_xy = (
            3 * x1 * y1 / 80
            - x1 * y2 / 80
            - x1 * y3 / 80
            - x1 * y4 / 80
            - x2 * y1 / 80
            + 3 * x2 * y2 / 80
            - x2 * y3 / 80
            - x2 * y4 / 80
            - x3 * y1 / 80
            - x3 * y2 / 80
            + 3 * x3 * y3 / 80
            - x3 * y4 / 80
            - x4 * y1 / 80
            - x4 * y2 / 80
            - x4 * y3 / 80
            + 3 * x4 * y4 / 80
        )
        bar_xz = (
            3 * x1 * z1 / 80
            - x1 * z2 / 80
            - x1 * z3 / 80
            - x1 * z4 / 80
            - x2 * z1 / 80
            + 3 * x2 * z2 / 80
            - x2 * z3 / 80
            - x2 * z4 / 80
            - x3 * z1 / 80
            - x3 * z2 / 80
            + 3 * x3 * z3 / 80
            - x3 * z4 / 80
            - x4 * z1 / 80
            - x4 * z2 / 80
            - x4 * z3 / 80
            + 3 * x4 * z4 / 80
        )
        bar_yz = (
            3 * y1 * z1 / 80
            - y1 * z2 / 80
            - y1 * z3 / 80
            - y1 * z4 / 80
            - y2 * z1 / 80
            + 3 * y2 * z2 / 80
            - y2 * z3 / 80
            - y2 * z4 / 80
            - y3 * z1 / 80
            - y3 * z2 / 80
            + 3 * y3 * z3 / 80
            - y3 * z4 / 80
            - y4 * z1 / 80
            - y4 * z2 / 80
            - y4 * z3 / 80
            + 3 * y4 * z4 / 80
        )

        for i in range(10):
            dx = x[i, 0] - xc[0]
            dy = x[i, 1] - xc[1]
            dz = x[i, 2] - xc[2]

            A[icell, i, 0] = 1
            A[icell, i, 1] = dx
            A[icell, i, 2] = dy
            A[icell, i, 3] = dz
            A[icell, i, 4] = dx ** 2 / 2 - bar_xx
            A[icell, i, 5] = dy ** 2 / 2 - bar_yy
            A[icell, i, 6] = dz ** 2 / 2 - bar_zz
            A[icell, i, 7] = dx * dy - bar_xy
            A[icell, i, 8] = dx * dz - bar_xz
            A[icell, i, 9] = dy * dz - bar_yz

    return A
