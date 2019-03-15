# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from ocellaris.utils import get_local


class FreeSurfaceLocatorImplDG0:
    def __init__(self, simulation, c, value):
        self.simulation = simulation
        self.field = c
        self.value = value
        self.cache = preprocess(simulation, c)

    def compute_crossing_points(self):
        return crossing_points_and_cells(self.simulation, self.field, self.value, self.cache)


def preprocess(simulation, field):
    """
    Store DOF and geometry info that will not change unless the mesh is
    updated (which is not handled in any way)
    """
    mesh = simulation.data['mesh']
    dofmap = field.function_space().dofmap()
    conFC = simulation.data['connectivity_FC']

    # Data to compute
    facet_data = {}
    cell_dofs = [None] * mesh.num_cells()
    is_ghost_cell = [False] * mesh.num_cells()

    cell_coords = numpy.zeros((2, 3), float)
    for facet in dolfin.facets(mesh, 'all'):
        fid = facet.index()
        cell_ids = conFC(fid)
        if len(cell_ids) != 2:
            continue

        # Get midpoint coordinates of the two connected cells and the DG0
        # dof needed to find the field value in each of the cells
        for i, cell_id in enumerate(cell_ids):
            cell = dolfin.Cell(mesh, cell_id)
            is_ghost_cell[cell_id] = cell.is_ghost()

            pt = cell.midpoint()
            cell_coords[i, 0] = pt.x()
            cell_coords[i, 1] = pt.y()
            cell_coords[i, 2] = pt.z()

            # Get the one and only DG0 dof for this cell
            dofs = dofmap.cell_dofs(cell_id)
            assert len(dofs) == 1
            cell_dofs[cell_id] = dofs[0]

        # Unit vector from cell 1 to cell 0
        uvec = cell_coords[0] - cell_coords[1]
        uvec /= (uvec[0] ** 2 + uvec[1] ** 2 + uvec[2] ** 2) ** 0.5

        # Store data for later use
        coords0 = numpy.array(cell_coords[0])
        coords1 = numpy.array(cell_coords[1])
        facet_data[fid] = cell_ids[0], cell_ids[1], coords0, coords1, uvec

    return facet_data, cell_dofs, is_ghost_cell


def crossing_points_and_cells(simulation, field, value, preprocessed):
    """
    Find cells that contain the value iso surface. This is done by
    connecting cell midpoints across facets and seing if the level set
    crosses this line. If it does, the point where it crosses, the cell
    containing the free surface crossing point and the vector from the low
    value cell to the high value cell is stored.

    The direction vector is made into a unit vector and then scaled by the
    value difference between the two cells. The vectors are computed in
    this way so that they can be averaged (for a cell with multiple
    crossing points) to get an approximate direction of increasing value
    (typically increasing density, meaning they point into the fluid in a
    water/air simulation). This is used such that the high value and the
    low value sides of the field can be approximately determined.

    The field is assumed to be piecewice constant (DG0)
    """
    facet_data, cell_dofs, is_ghost_cell = preprocessed
    all_values = get_local(field)

    # We define acronym LCCM: line connecting cell midpoints
    #   - We restrinct ourselves to LCCMs that cross only ONE facet
    #   - We number LLCMs by the index of the crossed facet

    # Find the crossing points where the contour crosses a LCCM
    crossing_points = {}
    for fid, fdata in facet_data.items():
        # Get preprocessed data
        cid0, cid1, coords0, coords1, uvec = fdata

        # Check for level crossing
        v0 = all_values[cell_dofs[cid0]]
        v1 = all_values[cell_dofs[cid1]]
        b1, b2 = v0 < value, v1 < value
        if (b1 and b2) or not (b1 or b2):
            # LCCM not crossed by contour
            continue

        # Find the location where the contour line crosses the LCCM
        fac = (v0 - value) / (v0 - v1)
        crossing_point = tuple((1 - fac) * coords0 + fac * coords1)

        # Scaled direction vector
        direction = uvec * (v0 - v1)

        # Find the cell containing the contour line
        surf_cid = cid0 if fac <= 0.5 else cid1
        is_ghost = is_ghost_cell[cid0] if fac <= 0.5 else is_ghost_cell[cid1]

        if not is_ghost:
            # Store the point and direction towards the high value cell
            crossing_points.setdefault(surf_cid, []).append((crossing_point, direction))

    return crossing_points
