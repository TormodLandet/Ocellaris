import numpy
import dolfin
from ocellaris.utils import get_local


class FreeSurfaceLocatorImplDG0:
    def __init__(self, simulation, c, value):
        self.simulation = simulation
        self.field = c
        self.value = value

    def compute_crossing_points(self):
        return crossing_points_and_cells(self.simulation, self.field, self.value)


def crossing_points_and_cells(simulation, field, value):
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
    mesh = simulation.data['mesh']
    all_values = get_local(field)
    dofmap = field.function_space().dofmap()

    # Mesh connectivities
    conFC = simulation.data['connectivity_FC']

    # We define acronym LCCM: line connecting cell midpoints
    #   - We restrinct ourselves to LCCMs that cross only ONE facet
    #   - We number LLCMs by the index of the crossed facet

    # Find the crossing points where the contour crosses a LCCM
    cell_coords = numpy.zeros((2, 3), float)
    cell_values = numpy.zeros(2, float)
    crossing_points = {}
    for facet in dolfin.facets(mesh):
        fid = facet.index()
        cell_ids = conFC(fid)
        if len(cell_ids) != 2:
            continue

        is_ghost_cell = [False, False]
        for i, cell_id in enumerate(cell_ids):
            cell = dolfin.Cell(mesh, cell_id)
            is_ghost_cell[i] = cell.is_ghost()

            # LCCM endpoint coordinates
            pt = cell.midpoint()
            cell_coords[i, 0] = pt.x()
            cell_coords[i, 1] = pt.y()
            cell_coords[i, 2] = pt.z()

            # LCCM endpoint values
            dofs = dofmap.cell_dofs(cell_id)
            assert len(dofs) == 1
            cell_values[i] = all_values[dofs[0]]

        b1, b2 = cell_values[0] < value, cell_values[1] < value
        if (b1 and b2) or not (b1 or b2):
            # LCCM not crossed by contour
            continue

        # Find the location where the contour line crosses the LCCM
        v1, v2 = cell_values
        fac = (v1 - value) / (v1 - v2)
        crossing_point = tuple((1 - fac) * cell_coords[0] + fac * cell_coords[1])

        # Scaled direction vector
        direction = cell_coords[0] - cell_coords[1]
        direction /= (direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2) ** 0.5
        direction *= v1 - v2

        # Find the cell containing the contour line
        surf_cid = cell_ids[0] if fac <= 0.5 else cell_ids[1]
        is_ghost = is_ghost_cell[0] if fac <= 0.5 else is_ghost_cell[1]

        if not is_ghost:
            # Store the point and direction towards the high value cell
            crossing_points.setdefault(surf_cid, []).append((crossing_point, direction))

    return crossing_points
