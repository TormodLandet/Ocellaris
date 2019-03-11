# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from ocellaris.utils import timeit, OcellarisError


@timeit
def get_dof_region_marks(simulation, V):
    """
    Given a function space, return a dictionary mapping dof number to a
    list of region number (indexed from 0, same as region list from the
    input file).

    Many dofs will not be included in the mapping since they are not
    inside a boundary region (not on a boundary facet). This property is
    used elsewhere to identify boundary dofs, in mark_cell_layers() and
    in SlopeLimiterBoundaryConditions
    """
    # This function only supports a small subset of function spaces
    family = V.ufl_element().family()
    assert family in ('Lagrange', 'Discontinuous Lagrange')

    # Get local indices for the facet dofs for each facet in the cell
    facet_dof_indices = get_facet_dof_indices(V)

    # Get dofs that share the same location (relevant for DG)
    same_loc_dofs = get_same_loc_dofs(V)

    # Loop over mesh and get dofs that are connected to boundary regions
    dm = V.dofmap()
    facet_marks = [int(m) for m in simulation.data['boundary_marker'].array()]
    mesh = simulation.data['mesh']
    dof_region_marks = {}
    for cell in dolfin.cells(mesh, 'all'):
        dofs = dm.cell_dofs(cell.index())

        for ifacet, facet in enumerate(dolfin.facets(cell)):
            # Get facet region marker
            mark = facet_marks[facet.index()] - 1

            # Skip non-boundary facets
            if mark == -1:
                continue

            facet_dofs = dofs[facet_dof_indices[ifacet]]
            for fdof in facet_dofs:
                dof_region_marks.setdefault(fdof, []).append(mark)

    # Treat all dofs in the same location in the same way
    for fdof, regions in list(dof_region_marks.items()):
        for dof2 in same_loc_dofs[fdof]:
            if dof2 not in dof_region_marks:
                dof_region_marks[dof2] = regions
                continue
            for mark in regions:
                if mark not in dof_region_marks[dof2]:
                    dof_region_marks[dof2].append(mark)

    # Order must be the same on all ranks independent of iteration order
    # of the dof_region_marks dictionary in the above loop
    for marks in dof_region_marks.values():
        marks.sort()

    return dof_region_marks


def get_facet_dof_indices(V):
    """
    Get local indices for the facet dofs for each facet in the cell
    """
    ndim = V.mesh().topology().dim()
    degree = V.ufl_element().degree()
    ndim_deg = (ndim, degree)

    if ndim_deg == (2, 1):
        # Linear triangle
        facet_dof_indices = numpy.zeros((3, 2), dtype=int)
        facet_dof_indices[0, :] = (1, 2)
        facet_dof_indices[1, :] = (0, 2)
        facet_dof_indices[2, :] = (0, 1)
    elif ndim_deg == (2, 2):
        # Quadratic triangle
        facet_dof_indices = numpy.zeros((3, 3), dtype=int)
        facet_dof_indices[0, :] = (1, 2, 3)
        facet_dof_indices[1, :] = (0, 2, 4)
        facet_dof_indices[2, :] = (0, 1, 5)
    elif ndim_deg == (3, 1):
        # Linear tetrahedron
        facet_dof_indices = numpy.zeros((4, 3), dtype=int)
        facet_dof_indices[0, :] = (1, 2, 3)
        facet_dof_indices[1, :] = (0, 2, 3)
        facet_dof_indices[2, :] = (0, 1, 3)
        facet_dof_indices[3, :] = (0, 1, 2)
    elif ndim_deg == (3, 2):
        # Quadratic tetrahedron
        facet_dof_indices = numpy.zeros((4, 6), dtype=int)
        facet_dof_indices[0, :] = (1, 2, 3, 4, 5, 6)
        facet_dof_indices[1, :] = (0, 2, 3, 4, 7, 8)
        facet_dof_indices[2, :] = (0, 1, 3, 5, 7, 9)
        facet_dof_indices[3, :] = (0, 1, 2, 6, 8, 9)
    else:
        raise OcellarisError(
            'Unsupported element ndim=%d degree=%d' % ndim_deg,
            'The boundary condition get_dof_region_marks '
            'code does not support this element',
        )

    return facet_dof_indices


def get_same_loc_dofs(V):
    """
    Return a dictionary mapping dof number to other dofs at the same
    location in space. V should obviously be a discontinuous space,
    otherwise there will not be multiple dofs in the same location
    """
    gdim = V.mesh().geometry().dim()
    dof_coordinates = V.tabulate_dof_coordinates().reshape((-1, gdim))

    # Map dof coordinate to dofs, this is for DG so multiple dofs
    # will share the same location
    coord_to_dofs = {}
    max_neighbours = 0
    for dof in range(len(dof_coordinates)):
        coord = tuple(round(x, 5) for x in dof_coordinates[dof])
        dofs = coord_to_dofs.setdefault(coord, [])
        dofs.append(dof)
        max_neighbours = max(max_neighbours, len(dofs) - 1)

    # Loop through dofs at same coordinate and map them to each other
    same_loc_dofs = {}
    for dofs in coord_to_dofs.values():
        for dof in dofs:
            same_loc_dofs[dof] = tuple(d for d in dofs if d != dof)

    return same_loc_dofs


def mark_cell_layers(
    simulation, V, layers=0, dof_region_marks=None, named_boundaries=None
):
    """
    Return all cells on the boundary and all connected cells in a given
    number of layers surrounding the boundary cells. Vertex neighbours
    are used to determine a cells neighbours.

    The initial list of cells is taken from a dictionary mapping dof to
    boundary region number. All cells containing these dofs are taken as
    the zeroth layer. If no such dictionary is provided then
    get_dof_region_marks(simulation, V) is used, hence the cells
    containing the boundary facing facet are used as the zeroth level.

    If named_boundaries is given then only cells on the given named
    boundaries are included in the zeroth layer. The special name 'all'
    will include all boundary regions.

    @return: set of the cell numbers of marked cells
    """
    if named_boundaries is None:
        named_boundaries = ['all']
    if not named_boundaries:
        return set()
    if dof_region_marks is None:
        dof_region_marks = get_dof_region_marks(simulation, V)

    # Get all regions with names
    all_regions = {region.name: region.index for region in simulation.data['boundary']}
    all_idxs = {region.index: region.name for region in simulation.data['boundary']}

    # Verify that the boundary names are unique
    assert len(all_regions) == len(simulation.data['boundary'])
    assert len(all_regions) == len(all_idxs)

    # Check that all named regions correspond to boundary regions
    for rname in named_boundaries:
        if rname != 'all' and rname not in all_regions:
            raise OcellarisError(
                'Unknown boundary region in input',
                '%r is not a boundary region' % rname,
            )

    # Names of all boundary regions
    if 'all' not in named_boundaries:
        to_keep = [mark for mark, name in all_idxs.items() if name in named_boundaries]

        # Remove marks to unwanted bounbdary regions
        drm = {
            dof: [mark for mark in dof_marks if mark in to_keep]
            for dof, dof_marks in dof_region_marks.items()
        }

        # Remove dofs wih empty region mark list
        dof_region_marks = {k: v for k, v in drm.items() if v}

    # Mark initial zeroth layer cells
    mesh = simulation.data['mesh']
    dm = V.dofmap()
    boundary_cells = set()
    for cell in dolfin.cells(mesh):
        dofs = dm.cell_dofs(cell.index())
        for dof in dofs:
            if dof in dof_region_marks:
                boundary_cells.add(cell.index())
                continue

    # Iteratively mark cells adjacent to boundary cells
    for _ in range(layers):
        boundary_cells_old = set(boundary_cells)
        for cell_index in boundary_cells_old:
            vertices = simulation.data['connectivity_CV'](cell_index)
            for vert_index in vertices:
                for nb in simulation.data['connectivity_VC'](vert_index):
                    boundary_cells.add(nb)

    return boundary_cells
