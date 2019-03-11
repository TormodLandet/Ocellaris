# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import ocellaris_interpolate, ocellaris_error
from . import register_multi_phase_model, MultiPhaseModel
from .vof import VOFMixin


# Default values, can be changed in the input file
CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION = False
MINIMUM_DIAMETER_OF_VERTEX_COLUMN = 1e-3


@register_multi_phase_model('HeightFunctionALE')
class HeightFunctionALE(VOFMixin, MultiPhaseModel):
    description = 'Separate phases with a height function and vertical mesh morphing'

    def __init__(self, simulation):
        """
        A simple height function multiphase model. The mesh is moved
        according to the vertical fluid velocity at the free surface
        after each time step.
        """
        self.simulation = simulation
        simulation.log.info('Creating height function mesh morphing multiphase model')

        # Define colour function
        V = simulation.data['Vc']
        c = dolfin.Function(V)
        simulation.data['c'] = c

        # Setup the mesh morpher
        simulation.mesh_morpher.setup()

        # Calculate mu from rho and nu (i.e mu is quadratic in c) or directly from c (linear in c)
        self.calculate_mu_directly_from_colour_function = simulation.input.get_value(
            'multiphase_solver/calculate_mu_directly_from_colour_function',
            CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION,
            'bool',
        )

        # The mean free surface height
        height_function_mean = simulation.input.get_value(
            'multiphase_solver/height_function_mean', required_type='float'
        )

        # The C++ code for the free surface height (final_height =
        # height_function_mean + height_function_cpp)
        height_function_cpp = simulation.input.get_value(
            'multiphase_solver/height_function_cpp', required_type='string'
        )
        assert 'x[1]' not in height_function_cpp and 'x[2]' not in height_function_cpp

        # What is the minimum width of a vertex column
        self.minimum_diameter_of_vertex_column = simulation.input.get_value(
            'multiphase_solver/minimum_diameter_of_vertex_column',
            MINIMUM_DIAMETER_OF_VERTEX_COLUMN,
            required_type='float',
        )

        # Morph the mesh to match the initial configuration of the free surface
        columns = initial_mesh_morphing(
            simulation,
            height_function_cpp,
            height_function_mean,
            eps=self.minimum_diameter_of_vertex_column,
        )
        self.vertex_columns = columns

        # Show some information about the columns
        simulation.log.info('    Created %d mesh columns' % len(columns))
        for col in columns:
            simulation.log.info('        %r' % col)

        # Get the physical properties
        self.rho0 = self.simulation.input.get_value(
            'physical_properties/rho0', required_type='float'
        )
        self.rho1 = self.simulation.input.get_value(
            'physical_properties/rho1', required_type='float'
        )
        self.nu0 = self.simulation.input.get_value(
            'physical_properties/nu0', required_type='float'
        )
        self.nu1 = self.simulation.input.get_value(
            'physical_properties/nu1', required_type='float'
        )

        # Update the rho and nu fields after each time step
        simulation.hooks.add_post_timestep_hook(
            self.update, 'HeightFunctionALE - update mesh'
        )

    def get_colour_function(self, k):
        """
        The colour function follows the cells and does not ever change
        """
        return self.simulation.data['c']

    def update(self):
        """
        Update the mesh position according to the calculated fluid velocities
        """
        timer = dolfin.Timer('Ocellaris HeightFunctionALE mesh update')
        sim = self.simulation

        # Get updated mesh velocity in each column
        # Use the fluid velocity at the free surface as the guiding velocity
        u_vertical = sim.data['u1']
        fs_vert_velocity = []
        for col in self.vertex_columns:
            dof = col.velocity_dof
            fs_vert_velocity.append(u_vertical.vector()[dof][0])

        # Move the mesh
        morph_mesh(self.simulation, self.vertex_columns, fs_vert_velocity)

        # Report properties of the colour field
        sum_c = dolfin.assemble(sim.data['c'] * dolfin.dx)
        sim.reporting.report_timestep_value('sum(c)', sum_c)

        timer.stop()


class MeshColumn(object):
    def __init__(self, start_pos, end_pos, fs_vert, fs_pos, fs_vel_dof):
        """
        Represent a column of mesh vertices. Each column has exactly one
        vertex on the free surface
        """
        # Info about the column
        self.start_pos = start_pos
        self.end_pos = end_pos

        # Info about the free surfaec vertex
        self.free_surface_vertex = fs_vert
        self.free_surface_pos = fs_pos
        self.velocity_dof = fs_vel_dof

        # Info about the rest of the vertices in the column
        self.vertices = []
        self.top = None
        self.bottom = None

    def __repr__(self, *args, **kwargs):
        return (
            '<MeshColumn: x in %r to %r y in %r to %r fs_vtx %r fs_pos %r vel_dof %r '
            'contains %d vertices'
        ) % (
            self.start_pos,
            self.end_pos,
            self.bottom,
            self.top,
            self.free_surface_vertex,
            self.free_surface_pos,
            self.velocity_dof,
            len(self.vertices),
        )


def initial_mesh_morphing(simulation, height_function_cpp, height_function_mean, eps):
    """
    This funcution does the following:

    - Identify the vertices closest to the free surface
    - Get the midpoints between such vertices. This defines a column
      with a startpos and enpos (in x) and a vertex index on the free
      surface.
    - Figure out how to move the free surface vertices to make them
      align with the free surface
    - Morph the mesh to fit the free surface with help from the defined
      columns
    - Return the column definitions for later use in mesh deformations
      as a list of :class:`MeshColumn` objects
    """
    D = simulation.ndim
    assert D == 2

    def get_height_func():
        height = ocellaris_interpolate(
            simulation, height_function_cpp, 'Height function', Vmesh
        )
        height.vector()[:] += height_function_mean
        return height

    mesh = simulation.data['mesh']
    Vmesh = simulation.data['Vmesh']
    Vu = simulation.data['Vu']
    Vc = simulation.data['Vc']
    con_CV = simulation.data['connectivity_CV']

    dofmap = Vmesh.dofmap()
    dofmap_fluid = Vu.dofmap()
    dofmap_colour = Vc.dofmap()
    colour_vec = simulation.data['c'].vector()
    height = get_height_func()

    # Find cells that contain the free surface
    vertices_to_move = set()
    vertex_coords = {}
    vertex_heights = {}
    vertex_dofs = {}
    vertex_fluid_vel_dofs = {}
    for cell in dolfin.cells(mesh):
        cid = cell.index()
        dofs = dofmap.cell_dofs(cid)
        dofs_fluid = dofmap_fluid.cell_dofs(cid)
        vertices = con_CV(cid)
        assert len(dofs) == len(vertices) == 3

        above = []
        below = []
        for i, vid in enumerate(vertices):
            if not vid in vertex_coords:
                vtx = dolfin.Vertex(mesh, vid)
                pos = [vtx.x(d) for d in range(D)]
                dof = dofs[i]
                h = height.vector()[dof][0]
                vertex_coords[vid] = pos
                vertex_heights[vid] = h
                vertex_dofs[vid] = dof
                vertex_fluid_vel_dofs[vid] = dofs_fluid[i]

            # h = vertex_heights[vid]
            h = height_function_mean
            above.append(vertex_coords[vid][1] >= h - 1e6)
            below.append(vertex_coords[vid][1] <= h + 1e6)

        if (all(above) and not any(below)) or (all(below) and not any(above)):
            # The free surface does not cut through this cell
            continue

        # The vertices to move to the free surface location
        # (We will end up moving only the closest vertex in each column)
        vertices_to_move.update(vertices)

    # Find x-positions of the vertices near the free surface
    xpos = [
        (vertex_coords[vid][0], vertex_coords[vid][1], vid) for vid in vertices_to_move
    ]
    assert len(xpos) > 1
    xpos.append((-1e10, None, None))
    xpos.append((1e10, None, None))
    xpos.sort()
    assert xpos[0][2] is None and xpos[-1][2] is None

    # Put the free surface vertices into columns. Make sure the columns
    # have a certain minimum width to be able to gather a few non-free
    # surface vertices in the column and to avoid moving more than one
    # vertex to the same location (in case of them being on the same
    # vertical line)
    columns = []
    in_col = []
    for i in range(1, len(xpos) - 1):
        x = xpos[i][0]
        in_col.append(xpos[i])

        # Get column start and end position
        if len(in_col) == 1:
            startpos = (
                (x + xpos[i - 1][0]) / 2
                if xpos[i - 1][2] is not None
                else x - eps / 100
            )
        endpos = (
            (x + xpos[i + 1][0]) / 2 if xpos[i + 1][2] is not None else x + eps / 100
        )

        # print 'x', x, 'xi-1', xpos[i-1][0], 'xi+1', xpos[i+1][0], 'y',
        # xpos[i][1], 'vid', xpos[i][2]
        if xpos[i + 1][0] - x < eps:
            continue
        # print 'NWCOL'

        # Get the closest y-coordinate (of free surface vertices) in the column
        min_diff_y = (None, None, None, 1e100)
        for p in in_col:
            x, y, vid = p
            dof = vertex_dofs[vid]
            # h = vertex_heights[vid]
            h = height_function_mean
            diff_y = abs(y - h)
            if diff_y < min_diff_y[3]:
                min_diff_y = (x, y, vid, diff_y)
        # print 'FS', min_diff_y

        # Create a MeshColumn object
        vid = min_diff_y[2]
        y_pos = min_diff_y[1]
        fluid_dof = vertex_fluid_vel_dofs[vid]
        col = MeshColumn(startpos, endpos, vid, y_pos, fluid_dof)
        columns.append(col)

        # Reset column
        in_col = []

    assert len(in_col) == 0

    # Put all vertices into columns
    for vid, coords in vertex_coords.items():
        for col in columns:
            if col.start_pos <= coords[0] < col.end_pos:
                col.vertices.append((vid, coords, vertex_dofs[vid]))
                break
        else:
            ocellaris_error(
                'HeightFunctionALE setup error',
                'Vertex at x=%f does not belong to any vertical free surface column'
                % coords[0],
            )

    # Find top and bottomn of each column
    for col in columns:
        miny, maxy = 1e10, -1e10
        for _, coords, _ in col.vertices:
            miny = min(coords[1], miny)
            maxy = max(coords[1], maxy)
        col.top = maxy
        col.bottom = miny

    ##############################################

    # Find the velocity with which to move the mesh to align it with the free surface
    vels = []
    for col in columns:
        vid = col.free_surface_vertex
        h = vertex_heights[vid]
        y = vertex_coords[vid][1]
        vels.append(h - y)

    # Move the mesh
    old_dt = simulation.dt
    simulation.dt = 1.0
    morph_mesh(simulation, columns, vels)
    simulation.dt = old_dt

    # Reset the mesh velocities
    for d in range(simulation.ndim):
        simulation.data['u_mesh%d' % d].vector().zero()

    ##############################################

    # Initialize the colour field. This will not change during the simulation.
    # Wet cells will remain wet and vice versa
    height = get_height_func()
    for cell in dolfin.cells(mesh):
        cid = cell.index()
        dofs = dofmap.cell_dofs(cid)
        h = sum(height.vector()[dof][0] for dof in dofs) / len(dofs)
        y = cell.midpoint().y()

        dofs_colour = dofmap_colour.cell_dofs(cid)
        assert len(dofs_colour) == 1
        dof_colour = dofs_colour[0]
        if y > h:
            colour_vec[dof_colour] = 0.0
        else:
            colour_vec[dof_colour] = 1.0

    return columns


def morph_mesh(simulation, columns, fs_vert_velocity):
    """
    Move the mesh. Each column is given a velocity in the
    vertical direction. All vertices in this column is moved
    according to the
    """
    assert len(columns) == len(fs_vert_velocity)
    dt = simulation.dt
    u_mesh = simulation.data['u_mesh1']

    for col, vel in zip(columns, fs_vert_velocity):
        y_fs = col.free_surface_pos
        for _vid, coords, dof in col.vertices:
            y_vtx = coords[1]
            if y_vtx > y_fs:
                fac = (col.top - y_vtx) / (col.top - y_fs)
            else:
                fac = (y_vtx - col.bottom) / (y_fs - col.bottom)
            u_mesh.vector()[dof] = vel * fac
        col.free_surface_pos = y_fs + vel * dt

        if col.free_surface_pos >= col.top:
            ocellaris_error(
                'HeightFunctionALE morphing error',
                'Free surface is above top of column in column %r' % col,
            )

    simulation.mesh_morpher.morph_mesh()
