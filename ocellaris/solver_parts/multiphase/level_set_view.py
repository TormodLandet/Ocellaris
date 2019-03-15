# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from collections import deque
import numpy
import dolfin
from ocellaris.probes.free_surface_locator import get_free_surface_locator
from ocellaris.utils import get_local


class LevelSetView:
    def __init__(self, simulation):
        """
        A LevelSetView is a view of another multi-phase model (typically
        VOF) as a level set function. This can be handy when distance from
        the free surface is needed in the code.

        This implementation of a level set function is not a stand-alone
        multi-phase model and cannot itself advect the free surface /
        density field.
        """
        self.simulation = simulation
        self.base = None
        self.base_type = 'unset'
        self.base_name = None
        self.iso_value = None
        self.name = None

        # Pieces of code that wants to know when we are updated
        self._callbacks = []

        # Create the level set function
        mesh = simulation.data['mesh']
        V = dolfin.FunctionSpace(mesh, 'CG', 1)
        self.level_set_function = dolfin.Function(V)
        self.cache = preprocess(simulation, self.level_set_function)

    def add_update_callback(self, cb):
        """
        Other functionality may depend on the level set view and want to
        be updated when we are updated
        """
        self._callbacks.append(cb)

    def set_density_field(self, c, name='c', value=0.5, update_hook='MultiPhaseModelUpdated'):
        """
        Set the generalised density field, either the vof function or
        the direct tensity field from a variable density simulation.
        This must be done before using the view
        """
        self.base = c
        self.base_type = 'vof'
        self.base_name = name
        self.iso_value = value

        # Store the level set function in an easy to access location and make
        # sure it is saved to any restart files that are written
        self.name = 'ls_%s_%s' % (name, float_to_ident(value))
        self.level_set_function.rename(self.name, self.name)
        self.simulation.data[self.name] = self.level_set_function

        # The locator helps find the position of the free surface
        self._locator = get_free_surface_locator(self.simulation, name, c, value)
        self._locator.add_update_hook(update_hook, self.update, 'Update LevelSetView')

        # Optionally add the level set view to plot output
        if self.simulation.input.get_value('multiphase_solver/plot_level_set_view', False, 'bool'):
            self.simulation.io.add_extra_output_function(self.level_set_function)

    def update(self):
        """
        The underlying data has changed, update the level set function
        """
        if self.base_type == 'vof':
            self._update_from_vof()
        else:
            raise NotImplementedError(
                'Cannot compute level set function ' 'from %r base field' % self.base_type
            )

        # Inform dependent functionality that we have updated
        for cb in self._callbacks:
            cb()

    def _update_from_vof(self):
        # This can be expensive, will involve recomputing the crossing
        # points if the density function has changed since the last access
        crossings = self._locator.crossing_points
        update_level_set_view(self.simulation, self.level_set_function, crossings, self.cache)


def float_to_ident(v):
    s = repr(v)
    for s0, s1 in [('.', '_'), ('+', 'p'), ('-', 'm')]:
        s = s.replace(s0, s1)
    return s


def preprocess(simulation, level_set_view):
    """
    Compute distance between dofs
    """
    V = level_set_view.function_space()
    mesh = V.mesh()
    dm = V.dofmap()

    # Get coordinates of both regular and ghost dofs
    dofs_x, Nlocal = all_dof_coordinates(V)

    cell_dofs = [None] * mesh.num_cells()
    dof_cells = [[] for _ in range(V.dim())]
    dof_dist = {}
    for cell in dolfin.cells(mesh, 'all'):
        cid = cell.index()
        dofs = dm.cell_dofs(cid)
        cell_dofs[cid] = list(dofs)
        for dof in dofs:
            dof_cells[dof].append(cid)

            # Store distance between the cell dofs
            for dof2 in dofs:
                if dof == dof2:
                    continue
                p1 = dofs_x[dof]
                p2 = dofs_x[dof2]
                vec = p1 - p2
                d = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5
                dof_dist[(dof, dof2)] = d
                dof_dist[(dof2, dof)] = d

    return dofs_x, dof_dist, cell_dofs, dof_cells


def update_level_set_view(simulation, level_set_view, crossings, cache):
    """
    Create a level set CG1 scalar function where the value is 0 at the
    given crossing point locations and approximately the distance to the
    nearest crossing point by following the edges of the mesh away from
    the crossing points and keeping track of the closest such point
    """
    dofs_x, dof_dist, cell_dofs, dof_cells = cache

    values = level_set_view.vector().get_local()
    values[:] = 1e100
    Nlocal = len(values)

    # For MPI ranks that contain the free surface we first fill out the
    # distance values starting at the free surface
    if crossings:
        # Mark distances in cells with a free surface
        for cid, cross in crossings.items():
            for crossing_point, direction in cross:
                for dof in cell_dofs[cid]:
                    if dof < Nlocal:
                        dpos = dofs_x[dof]
                        vec = dpos - crossing_point
                        dist = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5

                        # Make the level set function a signed distance function
                        # if vec.dot(direction) < 0:
                        #     dist = -dist

                        if abs(dist) < abs(values[dof]):
                            values[dof] = dist

        # Breadth-first search to populate all of the local domain
        queue = deque()
        for cid in crossings:
            queue.extend([d for d in cell_dofs[cid] if d < Nlocal])
        bfs(queue, values, cell_dofs, dof_cells, dof_dist)

    # Update ghost cell values
    level_set_view.vector().set_local(values)
    level_set_view.vector().apply('insert')

    # Breadth first search from ghost dofs
    values2 = get_local(level_set_view)
    queue = deque(range(Nlocal, len(values2)))
    bfs(queue, values2, cell_dofs, dof_cells, dof_dist)
    level_set_view.vector().set_local(values2[:Nlocal])
    level_set_view.vector().apply('insert')


def bfs(queue, values, cell_dofs, dof_cells, dof_dist):
    """
    Breadth-first search to populate all of the local domain, starting
    from a deque "queue" which contains dofs indices where the search
    should start
    """
    checks = 0
    Nval = len(values)

    while queue:
        dof = queue.popleft()
        checks += 1
        dval = values[dof]

        for cid in dof_cells[dof]:
            for dof2 in cell_dofs[cid]:
                if dof2 == dof or dof2 >= Nval:
                    continue

                # Compute the additional distance to dof2
                dx = dof_dist[(dof, dof2)]
                dval2 = values[dof2]

                # Used when the level set function is the signed distance.
                # Has some problems with inaccurate surface normals
                # s = 1 if dval >= 0 else -1
                # distance_to_crossing = dval + s * dx

                # Unsigned distance function
                distance_to_crossing = dval + dx

                # Update if we have found a shorter path to a crossing
                if abs(distance_to_crossing) < abs(dval2):
                    values[dof2] = distance_to_crossing

                    # Update neighbours since we changed
                    queue.append(dof2)

    return checks


def all_dof_coordinates(V):
    """
    Return coordinates of all dofs, similar to tabulate_dof_coordinates,
    but also includes the ghosts dofs. This always returns 3-vectors for
    the coordinates, even if the mesh is 2D or 1D.

    This has all communication going via the root. Consider rewriting this
    by using DofMap.off_process_owner to avoid the root botleneck for
    large MPI runs
    """
    mesh = V.mesh()

    # Get coordinates of regular dofs
    gdim = mesh.geometry().dim()
    dofs_x = V.tabulate_dof_coordinates().reshape((-1, gdim))
    global_dofs = V.dofmap().tabulate_local_to_global_dofs()
    Nlocal = len(dofs_x)

    # Make sure dof positions are 3-vectors
    if gdim != 3:
        tmp = numpy.zeros((dofs_x.shape[0], 3), float)
        tmp[:, :gdim] = dofs_x
        dofs_x = tmp

    # Do nothing in serial
    comm = mesh.mpi_comm()
    if comm.size == 1:
        return dofs_x, Nlocal

    # Missing and present information
    missing_dofs = global_dofs[Nlocal:]
    present_dofs = global_dofs[:Nlocal]

    # Send list of global dofs which we want to know the coordinates for
    # along with the global dofs where we know the coordinates
    alldata = comm.gather((missing_dofs, present_dofs, dofs_x))

    # Global dof coordinate map, may need to be changed for comm.size >> 1
    if comm.rank == 0:
        # Build global dof coordinate mapping using each ranks info
        global_dofs_x = {}
        for _rank_missing_dofs, rank_present_dofs, rank_dofs_x in alldata:
            for gdof, coord in zip(rank_present_dofs, rank_dofs_x):
                global_dofs_x[gdof] = coord

        # Find the missing local data in the global mapping for each rank
        ret_data = []
        for rank_missing_dofs, _rank_present_dofs, _rank_dofs_x in alldata:
            ret = []
            ret_data.append(ret)
            for gdof in rank_missing_dofs:
                ret.append(global_dofs_x[gdof])
    else:
        ret_data = None

    # Get back the missing values
    ghost_coords = comm.scatter(ret_data)

    # Create dof coordinates with ghost positions
    dofs_x_all = numpy.zeros((len(global_dofs), 3), float)
    dofs_x_all[:Nlocal] = dofs_x
    for i, coord in enumerate(ghost_coords):
        dofs_x_all[Nlocal + i] = coord

    return dofs_x_all, Nlocal

