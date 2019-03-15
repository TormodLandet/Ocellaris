# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
import numpy
from matplotlib import pyplot
from ocellaris.utils import gather_lines_on_root, timeit, ocellaris_error
from . import Probe, register_probe


INCLUDE_BOUNDARY = False
WRITE_INTERVAL = 1
SHOW_INTERVAL = 0


@register_probe('IsoSurface')
class IsoSurface(Probe):
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.cells_with_surface = None
        inp = probe_input

        assert (
            self.simulation.ndim == 2
        ), 'IsoSurface only implemented in 2D (contour line)'
        assert not self.simulation.mesh_morpher.active, 'IsoSurface does not support ALE yet'

        # Read input
        self.name = inp.get_value('name', required_type='string')
        self.field_name = inp.get_value('field', required_type='string')
        self.value = inp.get_value('value', required_type='float')
        self.custom_hook_point = inp.get_value(
            'custom_hook', None, required_type='string'
        )
        self.field = simulation.data[self.field_name]
        self.include_boundary = inp.get_value(
            'include_boundary', INCLUDE_BOUNDARY, 'bool'
        )

        # Should we write the data to a file
        prefix = simulation.input.get_value('output/prefix', None, 'string')
        file_name = inp.get_value('file_name', None, 'string')
        self.write_file = file_name is not None
        if self.write_file:
            if prefix is not None:
                self.file_name = prefix + file_name
            else:
                self.file_name = file_name
            self.write_interval = inp.get_value('write_interval', WRITE_INTERVAL, 'int')

        # Should we pop up a matplotlib window when running?
        self.show_interval = inp.get_value('show_interval', SHOW_INTERVAL, 'int')
        self.show = self.show_interval != 0 and simulation.rank == 0
        self.xlim = inp.get_value('xlim', (None, None), 'list(float)')
        self.ylim = inp.get_value('ylim', (None, None), 'list(float)')
        if not len(self.xlim) == 2:
            ocellaris_error(
                'Plot xlim must be two numbers',
                'IsoSurface probe "%s" contains invalid xlim specification' % self.name,
            )
        if not len(self.ylim) == 2:
            ocellaris_error(
                'Plot ylim must be two numbers',
                'IsoSurface probe "%s" contains invalid ylim specification' % self.name,
            )

        if self.write_file and simulation.rank == 0:
            self.output_file = open(self.file_name, 'wt')
            self.output_file.write(
                '# Ocellaris iso surface of the %s field\n' % self.field_name
            )
            self.output_file.write('# value = %15.5e\n' % self.value)
            self.output_file.write('# dim = %d\n' % self.simulation.ndim)

        if self.show and simulation.rank == 0:
            pyplot.ion()
            self.fig = pyplot.figure()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_title('Iso surface %s' % self.name)

        # The class that finds the contour line
        if self.field_name == 'height_function':
            self.locator = HeightFunctionLocator(simulation)
        else:
            self.locator = IsoSurfaceLocator(simulation, self.field.function_space())

        if self.custom_hook_point is not None:
            simulation.hooks.add_custom_hook(
                self.custom_hook_point, self.run, 'Probe "%s"' % self.name
            )
        else:
            self.end_of_timestep = self.run

        # Flush output file
        self.simulation.hooks.add_custom_hook('flush', self.flush, 'Flush log file')

    def run(self, force_active=False):
        """
        Find and output the line probe
        """
        it = self.simulation.timestep

        # Should we update the plot?
        update_plot = False
        if self.show and (it == 1 or it % self.show_interval == 0):
            update_plot = True

        # Should we update the file?
        update_file = False
        if self.write_file and (it == 1 or it % self.write_interval == 0):
            update_file = True

        # Do not do any postprocessing for non-requested time steps
        if not (update_file or update_plot or force_active):
            return

        # Get the iso surfaces
        surfaces, cells = self.locator.get_iso_surface(self.field, self.value)
        self.cells_with_surface = cells

        # Get the boundary surfaces for cells with values above the given value
        if self.include_boundary:
            surfaces += get_boundary_surface(self.simulation, self.field, self.value)

        # Create lines (this assumes 2D and throws away the z-component)
        lines = []
        for surface in surfaces:
            x = numpy.array([pos[0] for pos in surface], float)
            y = numpy.array([pos[1] for pos in surface], float)
            lines.append((x, y))

        # Communicate lines to the root process in case we are running in parallel
        gather_lines_on_root(lines)

        if update_file and self.simulation.rank == 0:
            self.output_file.write(
                'Time %10.5f nsurf %d\n' % (self.simulation.time, len(lines))
            )
            for x, y in lines:
                self.output_file.write(' '.join('%10.5f' % v for v in x) + '\n')
                self.output_file.write(' '.join('%10.5f' % v for v in y) + '\n')
                self.output_file.write(' '.join('%10.5f' % 0 for v in x) + '\n')

        if update_plot and self.simulation.rank == 0:
            with dolfin.Timer('Ocellaris plot iso surface'):
                self.ax.clear()
                for x, y in lines:
                    self.ax.plot(x, y)
                self.ax.set_xlabel('x')
                self.ax.set_ylabel('y')
                self.ax.relim()
                self.ax.autoscale_view()
                if self.xlim != (None, None):
                    self.ax.set_xlim(*self.xlim)
                if self.ylim != (None, None):
                    self.ax.set_ylim(*self.ylim)
                self.fig.canvas.draw()
                # self.fig.canvas.flush_events()

        # Return value only used in unit testing
        return lines

    def flush(self):
        if (
            self.write_file
            and self.simulation.rank == 0
            and not self.output_file.closed
        ):
            self.output_file.flush()

    def end_of_simulation(self):
        """
        The simulation is done. Close the output file
        """
        if self.write_file and self.simulation.rank == 0:
            self.output_file.close()


class IsoSurfaceLocatorCache:
    pass


class IsoSurfaceLocator:
    def __init__(self, simulation, V):
        self.simulation = simulation
        self.degree = V.ufl_element().degree()

        self.cache = IsoSurfaceLocatorCache()
        if self.degree in (1, 2):
            return prepare_DG1_DG2(self.cache, simulation, V)

    @timeit
    def get_iso_surface(self, field, value):
        """
        Find the iso-surfaces (contour lines) of the
        given field with the given scalar value
        """
        sim = self.simulation
        if self.degree == 0:
            return get_iso_surfaces_picewice_constants(sim, field, value)
        elif self.degree in (1, 2):
            return get_iso_surface_DG1_DG2(sim, self.cache, field, value)
        else:
            return get_iso_surfaces(sim, field, value)


class HeightFunctionLocator:
    def __init__(self, simulation):
        self.simulation = simulation

    @timeit
    def get_iso_surface(self, field, value):
        """
        Find the iso-surfaces (contour lines) of the
        given field with the given scalar value
        """
        sim = self.simulation
        xpos = sim.data['height_function_x']
        ypos = (
            sim.data['height_function'].vector().get_local()
        )  # only needs the first elements
        line = [(x, h) for x, h in zip(xpos, ypos)]
        return [line], None


def get_iso_surfaces(simulation, field, value):
    """
    Slow fallback version that uses vertex values only
    """
    assert simulation.ndim == 2
    mesh = simulation.data['mesh']
    all_values = field.compute_vertex_values()

    # We will collect the cells containing the iso surface
    cells_with_surface = numpy.zeros(mesh.num_cells(), bool)
    connectivity_FC = simulation.data['connectivity_FC']

    # Find the crossing points where the contour crosses a facet
    crossing_points = {}
    for facet in dolfin.facets(mesh):
        fid = facet.index()

        # Get connected vertices and the field values there
        vertex_coords = []
        vertex_values = []
        for vertex in dolfin.vertices(facet):
            pt = vertex.point()
            vertex_coords.append((pt.x(), pt.y(), pt.z()))
            vertex_values.append(all_values[vertex.index()])
        assert len(vertex_coords) == 2

        # Check for iso surface crossing
        b1, b2 = vertex_values[0] < value, vertex_values[1] < value
        if (b1 and b2) or not (b1 or b2):
            # Facet not crossed by contour
            continue

        # Find the location where the contour line crosses the facet
        v1, v2 = vertex_values
        fac = (v1 - value) / (v1 - v2)
        x = (1 - fac) * vertex_coords[0][0] + fac * vertex_coords[1][0]
        y = (1 - fac) * vertex_coords[0][1] + fac * vertex_coords[1][1]
        z = (1 - fac) * vertex_coords[0][2] + fac * vertex_coords[1][2]
        crossing_points[fid] = (x, y, z)

        # Find the cells connected to this facet
        for cid in connectivity_FC(fid):
            cells_with_surface[cid] = True

    # Get facet-facet connectivity via cells
    conFC = simulation.data['connectivity_FC']
    conCF = simulation.data['connectivity_CF']

    # Find facet to facet connections
    connections = {}
    for facet_id in crossing_points:
        connections[facet_id] = []
        for cell_id in conFC(facet_id):
            for facet_neighbour_id in conCF(cell_id):
                if (
                    facet_neighbour_id != facet_id
                    and facet_neighbour_id in crossing_points
                ):
                    connections[facet_id].append(facet_neighbour_id)

    # Make continous contour lines
    # Find end points of contour lines and start with these
    end_points = [
        facet_id for facet_id, neighbours in connections.items() if len(neighbours) == 1
    ]
    contours_from_endpoints = contour_lines_from_endpoints(
        end_points, crossing_points, connections
    )

    # Include crossing points without neighbours or joined circles without end points
    other_points = crossing_points.keys()
    contours_from_singles_and_loops = contour_lines_from_endpoints(
        other_points, crossing_points, connections
    )

    assert len(crossing_points) == 0
    return contours_from_endpoints + contours_from_singles_and_loops, cells_with_surface


def get_iso_surfaces_picewice_constants(simulation, field, value):
    """
    Find the iso-surfaces (contour lines) of the
    given field with the given scalar value

    The field is assumed to be piecewice constant (DG0)
    """
    assert simulation.ndim == 2
    mesh = simulation.data['mesh']
    all_values = field.vector().get_local()
    dofmap = field.function_space().dofmap()

    # Mesh connectivities
    conFC = simulation.data['connectivity_FC']
    conVF = simulation.data['connectivity_VF']
    conFV = simulation.data['connectivity_FV']

    # We will collect the cells containing the iso surface
    cells_with_surface = numpy.zeros(mesh.num_cells(), bool)

    # We define acronym LCCM: line connecting cell midpoints
    #   - We restrinct ourselves to LCCMs that cross only ONE facet
    #   - We number LLCMs by the index of the crossed facet

    # Find the crossing points where the contour crosses a LCCM
    vertex_coords = numpy.zeros((2, 3), float)
    vertex_values = numpy.zeros(2, float)
    crossing_points = {}
    for facet in dolfin.facets(mesh):
        fid = facet.index()
        cell_ids = conFC(fid)
        if len(cell_ids) != 2:
            continue

        has_ghost_cell = False
        for i, cell_id in enumerate(cell_ids):
            cell = dolfin.Cell(mesh, cell_id)
            if cell.is_ghost():
                has_ghost_cell = True
                break

            # LCCM endpoint coordinates
            pt = cell.midpoint()
            vertex_coords[i, 0] = pt.x()
            vertex_coords[i, 1] = pt.y()
            vertex_coords[i, 2] = pt.z()

            # LCCM endpoint values
            dofs = dofmap.cell_dofs(cell_id)
            assert len(dofs) == 1
            vertex_values[i] = all_values[dofs[0]]

        if has_ghost_cell:
            continue

        b1, b2 = vertex_values[0] < value, vertex_values[1] < value
        if (b1 and b2) or not (b1 or b2):
            # LCCM not crossed by contour
            continue

        # Find the location where the contour line crosses the LCCM
        v1, v2 = vertex_values
        fac = (v1 - value) / (v1 - v2)
        crossing_points[fid] = tuple(
            (1 - fac) * vertex_coords[0] + fac * vertex_coords[1]
        )

        # Find the cell containing the contour line
        surf_cid = cell_ids[0] if fac <= 0.5 else cell_ids[1]
        cells_with_surface[surf_cid] = True

    # Find facet to facet connections
    connections = {}
    for facet_id in crossing_points:
        connections[facet_id] = []
        for vertex_id in conFV(facet_id):
            for facet_neighbour_id in conVF(vertex_id):
                if (
                    facet_neighbour_id != facet_id
                    and facet_neighbour_id in crossing_points
                ):
                    connections[facet_id].append(facet_neighbour_id)

    # Make continous contour lines
    # Find end points of contour lines and start with these
    end_points = [
        facet_id for facet_id, neighbours in connections.items() if len(neighbours) == 1
    ]
    contours_from_endpoints = contour_lines_from_endpoints(
        end_points, crossing_points, connections
    )

    # Include crossing points without neighbours or joined circles without end points
    other_points = crossing_points.keys()
    contours_from_singles_and_loops = contour_lines_from_endpoints(
        other_points, crossing_points, connections
    )

    assert len(crossing_points) == 0
    return contours_from_endpoints + contours_from_singles_and_loops, cells_with_surface


def prepare_DG1_DG2(cache, simulation, V):
    """
    Prepare to find iso surfaces of the given field. Caches geometry and
    topology data of the mesh and must be rerun if this data changes!

    This is rather slow, but it runs only once and the results are cached
    """
    mesh = simulation.data['mesh']
    gdim = mesh.geometry().dim()
    degree = V.ufl_element().degree()
    dofmap = V.dofmap()
    dofs_x = V.tabulate_dof_coordinates().reshape((-1, gdim))
    N = dofs_x.shape[0]
    assert degree in (1, 2)
    assert gdim == 2

    # Map location to dof
    x_to_dof = {}
    for dof, pos in enumerate(dofs_x):
        x_to_dof.setdefault(tuple(pos), []).append(dof)

    # Create mapping from dof to other dofs at the same coordinate
    same_loc = [[] for _ in range(N)]
    for dofs in x_to_dof.values():
        for d0 in dofs:
            for d1 in dofs:
                if d0 != d1:
                    same_loc[d0].append(d1)

    # Find the immediate neighbours for all dofs where
    # immediate neighbour is defined as
    # either 1) sits in the same location (but different cell)
    # or     2) sits next to each other on the same facet
    immediate_neighbours = [None] * N
    connected_cells = [None] * N
    for cell in dolfin.cells(mesh):
        cid = cell.index()
        dofs = dofmap.cell_dofs(cid)
        if degree == 1:
            for dof in dofs:
                # Immediate neighbours
                nbs = list(same_loc[dof])
                immediate_neighbours[dof] = nbs
                nbs.extend(d for d in dofs if dof != d)

                # The first connected cell is the cell owning the dof
                connected_cells[dof] = [cid]
        else:
            for i, dof in enumerate(dofs):
                # Immediate neighbours
                nbs = list(same_loc[dof])
                immediate_neighbours[dof] = nbs
                if i == 0:
                    nbs.extend((dofs[4], dofs[5]))
                elif i == 1:
                    nbs.extend((dofs[3], dofs[5]))
                elif i == 2:
                    nbs.extend((dofs[3], dofs[4]))
                elif i == 3:
                    nbs.extend((dofs[1], dofs[2]))
                elif i == 4:
                    nbs.extend((dofs[0], dofs[2]))
                elif i == 5:
                    nbs.extend((dofs[0], dofs[1]))

                # The first connected cell is the cell owning the dof
                connected_cells[dof] = [cid]

    # Extend list of connected cells
    for dof, pos in enumerate(dofs_x):
        p = tuple(pos)
        for nb in x_to_dof[p]:
            if nb != dof:
                # Get the first connected cell (the cell containing the dof)
                nb_cid = connected_cells[nb][0]
                connected_cells[dof].append(nb_cid)

    # Find the extended neighbour dofs of all dofs. An extended neighbbour is
    # a dof that can be the next point on a contour line. The line between a
    # dof and its extended neighbour can hence not cut through a facet, but it
    # can be parallell/incident with a facet
    extended_neighbours = [None] * N
    for dof, cell_ids in enumerate(connected_cells):
        extended_neighbours[dof] = enbs = []
        for cid in cell_ids:
            # Add dofs in conneced cell as neighbours
            for d in dofmap.cell_dofs(cid):
                if d != dof and d not in enbs:
                    enbs.append(d)
                # Add other dofs at the same location
                for d2 in same_loc[d]:
                    if d2 != dof and d2 not in enbs:
                        enbs.append(d2)

    # Sort extended neighbours by distance
    for dof in range(N):
        enbs = extended_neighbours[dof]
        p = dofs_x[dof]
        tmp = []
        for n in enbs:
            pn = dofs_x[n]
            d = p - pn
            tmp.append((d.dot(d), n))
        tmp.sort(reverse=True)
        extended_neighbours[dof] = [n for _dist, n in tmp]

    cache.N = N
    cache.degree = degree
    cache.dofs_x = dofs_x
    cache.x_to_dof = x_to_dof
    cache.immediate_neighbours = immediate_neighbours
    cache.extended_neighbours = extended_neighbours
    cache.connected_cells = connected_cells


def get_iso_surface_DG1_DG2(simulation, cache, field, value):
    """
    Find the iso-surfaces (contour lines) of the given field with the
    given scalar value. The field is assumed to be linear or quadratic

    We assume that the field is discontinuous at internal facets. This
    means that the iso surface could be incident with a facet since
    the scalar field is dual valued there
    """
    all_values = field.vector().get_local()
    assert simulation.ndim == 2
    assert len(all_values) == cache.N

    # We will collect the cells containing the iso surface
    mesh = simulation.data['mesh']
    cells_with_surface = numpy.zeros(mesh.num_cells(), bool)

    # Find where the iso surface crosses between two dofs
    # Could be right at the same location, in the jump
    # between two cells
    crossing_points = {}
    for dof, nbs in enumerate(cache.immediate_neighbours):
        p0 = tuple(cache.dofs_x[dof])
        v0 = all_values[dof]
        for n in nbs:
            v1 = all_values[n]
            b0 = v0 <= value and v1 >= value
            b1 = v0 >= value and v1 <= value
            if not (b0 or b1):
                continue
            p1 = tuple(cache.dofs_x[n])

            # Check for surface at jump
            if p0 == p1:
                if dof == cache.x_to_dof[p0][0]:
                    crossing_points[dof] = p0

                    # Mark cell as surface cell
                    cid = cache.connected_cells[dof][0]
                    cells_with_surface[cid] = True

            # Check for surface crossing a facet
            elif b0:
                fac = (v0 - value) / (v0 - v1)
                cp = ((1 - fac) * p0[0] + fac * p1[0], (1 - fac) * p0[1] + fac * p1[1])

                # Select the dof with the fewest connections to avoid
                # the iso surface line crossing a facet (the dof with
                # the most extended neighbours will be the corner dof)
                num_enbs0 = len(cache.extended_neighbours[dof])
                num_enbs1 = len(cache.extended_neighbours[n])
                if num_enbs0 < num_enbs1:
                    crossing_points[dof] = cp
                else:
                    crossing_points[n] = cp

                # Mark cell as surface cell
                cid = cache.connected_cells[dof][0]
                cells_with_surface[cid] = True

    # Get connections between crossing points using the extended
    # dof neighbourhood to look for possible connections
    connections = {}
    for dof, pos in crossing_points.items():
        tmp = []
        for n in cache.extended_neighbours[dof]:
            if n not in crossing_points:
                continue
            p2 = crossing_points[n]
            d = (pos[0] - p2[0]) ** 2 + (pos[1] - p2[1]) ** 2
            tmp.append((d, n))
        tmp.sort(reverse=True)
        connections[dof] = [n for _, n in tmp]

    # Make continous contour lines
    possible_starting_points = crossing_points.keys()
    contours = contour_lines_from_endpoints(
        possible_starting_points,
        crossing_points,
        connections,
        min_length=3 * cache.degree,
        backtrack_from_end=True,
        extend_endpoints=False,
    )

    return contours, cells_with_surface


@timeit
def get_boundary_surface(simulation, field, value):
    """
    Find the boundary surface consisting of facets with
    scalar field values greater than the given iso value
    """
    assert simulation.ndim == 2

    mesh = simulation.data['mesh']
    all_values = field.compute_vertex_values()

    connectivity_FC = simulation.data['connectivity_FC']

    # Find the crossing points where the contour crosses a facet
    connections = {}
    for facet in dolfin.facets(mesh):
        fid = facet.index()

        # Skip facets that are not on the boundary
        connected_cells = connectivity_FC(fid)
        if not len(connected_cells) == 1:
            continue

        # Get connected vertices and the field values there
        vertex_coords = []
        vertex_values = []
        for vertex in dolfin.vertices(facet):
            pt = vertex.point()
            vertex_coords.append((pt.x(), pt.y(), pt.z()))
            vertex_values.append(all_values[vertex.index()])
        assert len(vertex_coords) == 2

        # Check if all values are above the iso value
        if vertex_values[0] < value or vertex_values[1] < value:
            continue

        connections.setdefault(vertex_coords[0], []).append(vertex_coords[1])
        connections.setdefault(vertex_coords[1], []).append(vertex_coords[0])

    # Map coord to coord, just to be able to use the generic functionality in
    # contour_lines_from_endpoints which works on facet_id <-> coord mappings
    available_coords = {vc: vc for vc in connections}

    # Make continous contour lines
    # Find end points of contour lines and start with these
    end_points = [vc for vc, neighbours in connections.items() if len(neighbours) < 2]
    contours_from_endpoints = contour_lines_from_endpoints(
        end_points, available_coords, connections
    )

    # Include crossing points without neighbours or joined circles without end points
    other_points = available_coords.keys()
    contours_from_singles_and_loops = contour_lines_from_endpoints(
        other_points, available_coords, connections
    )

    assert len(available_coords) == 0
    return contours_from_endpoints + contours_from_singles_and_loops


def contour_lines_from_endpoints(
    endpoints,
    crossing_points,
    connections,
    min_length=3,
    backtrack_from_end=False,
    extend_endpoints=True,
):
    """
    Follow contour lines and create contours

    - endpoints: an iterable of crossing point IDs (could be facet ids)
    - crossing_points: a mapping from crossing point ID to position coordinates
    - connections: a mapping from ID to IDs
    """
    contours = []
    endpoint_ids = list(endpoints)

    while endpoint_ids:
        end_id = endpoint_ids.pop()
        if not end_id in crossing_points:
            # This has been taken by the other end
            continue

        # Make a new contour line
        contour = [crossing_points.pop(end_id)]

        # Loop over neighbours to the end of the contour
        queue = list(connections[end_id])
        prev = end_id
        has_backtracked = False
        while queue:
            nb_id = queue.pop()

            # Is this neighbour a possible next part of the contour line?
            if nb_id in crossing_points:
                # Found connection to use as next in contour
                cpoint = crossing_points.pop(nb_id)
                if cpoint != contour[-1]:
                    contour.append(cpoint)

                # Unused connections may be new end points (more than two connections)
                if extend_endpoints:
                    for candidate in queue:
                        if candidate in crossing_points:
                            endpoint_ids.append(candidate)

                queue = list(connections[nb_id])
                prev = nb_id

            # Is this the end of a loop?
            if (
                nb_id == end_id
                and prev in connections[end_id]
                and len(contour) > min_length - 1
                and not has_backtracked
            ):
                contour.append(contour[0])
                break

            # Backtrack from endpoint in case it was not a true endpoint
            if backtrack_from_end and not queue and not has_backtracked:
                contour = contour[::-1]
                queue = list(connections[end_id])
                has_backtracked = True

        contours.append(contour)

    return contours
