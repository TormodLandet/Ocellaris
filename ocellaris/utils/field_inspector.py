# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
import numpy


class FieldInspector(object):
    def __init__(
        self, xdmf_file_name, V, startpos=(-1e100, -1e100), endpos=(1e100, 1e100)
    ):
        """
        This class is meant for producing visualisations of DG fields,
        most notably DG2 fields. We subdevide each cell to be able to
        properly display higher order fields, and we also make each cell
        disconnected from the neighbours such that discontinuities are
        apparent

        If you have more than a few cells in your mesh you should specify
        startpos and endpos to only include cells in that rectangle,
        otherwise this class will be quite slow.
        """
        self.xdmf_file_name = xdmf_file_name
        self.function_space = V
        self.startpos = startpos
        self.endpos = endpos

        self.mesh = V.mesh()
        self.family = V.ufl_element().family()
        self.degree = V.ufl_element().degree()

        # This will not work in parallel
        assert dolfin.MPI.size(dolfin.MPI.comm_world) == 1

        # Only tested for DG2 for now
        assert self.mesh.geometry().dim() == 2
        assert self.mesh.topology().dim() == 2
        assert self.family == 'Discontinuous Lagrange'
        assert self.degree == 2

        self.selected_cells = get_selected_cells(self.mesh, startpos, endpos)
        self.subdivided_mesh, self.parent_cell = create_subdivided_mesh(
            self.mesh, self.selected_cells, self.degree
        )
        self.subdivided_function_space = dolfin.FunctionSpace(
            self.subdivided_mesh, 'CG', 1
        )
        self.subdivided_function = dolfin.Function(self.subdivided_function_space)

        self.xdmf = dolfin.XDMFFile(dolfin.MPI.comm_world, xdmf_file_name)

    def write_function(self, func, t):
        """
        Write a function to XDMF file. The function is from the original
        function space, but will be written in the subdivided function space
        with truely discontinuous elements
        """
        f = self.subdivided_function
        V = self.subdivided_function_space
        dm = V.dofmap()
        gdim = self.subdivided_mesh.geometry().dim()
        dof_coords = V.tabulate_dof_coordinates().reshape((-1, gdim))
        all_values = f.vector().get_local()
        f.rename(func.name(), func.name())

        vals = numpy.zeros(1, float)
        crds = numpy.zeros(2, float)

        # Evaluate the function at all subcell dof locations
        for subcell in dolfin.cells(self.subdivided_mesh):
            subcell_id = subcell.index()
            for dof in dm.cell_dofs(subcell_id):
                crds[:] = dof_coords[dof]
                parent_cell = dolfin.Cell(self.mesh, self.parent_cell[subcell_id])
                func.eval_cell(vals, crds, parent_cell)
                all_values[dof] = vals[0]

        f.vector().set_local(all_values)
        f.vector().apply('insert')
        self.xdmf.write(f, float(t))


def get_selected_cells(mesh, startpos, endpos):
    """
    Return a list of cells contained in the startpos-endpos rectangle
    """
    xstart, ystart = startpos
    xend, yend = endpos

    selected_cells = set()
    vertex_coords = mesh.coordinates()
    for cell in dolfin.cells(mesh):
        cell_vertices = cell.entities(0)
        for vid in cell_vertices:
            x, y = vertex_coords[vid]
            if xstart <= x <= xend and ystart <= y <= yend:
                selected_cells.add(cell.index())
                break

    return selected_cells


def create_subdivided_mesh(mesh, selected_cells, degree):
    """
    Create a subdivided mesh for the given cells that can
    represent a function of the given degree
    """
    num_div = degree + 2

    # Subdivide each of the selected cells
    vertex_coords = mesh.coordinates()
    subvertices = []
    subcells = []
    parent_cell = {}
    for cell in dolfin.cells(mesh):
        cid = cell.index()
        if cid not in selected_cells:
            continue

        # Create subdivided cells
        vid0, vid1, vid2 = cell.entities(0)
        x0 = tuple(vertex_coords[vid0, :])
        x1 = tuple(vertex_coords[vid1, :])
        x2 = tuple(vertex_coords[vid2, :])
        subdiv_cells = subdivide_cell(num_div, x0, x1, x2)

        # Number the new vertices and subcells. Connectivity to other
        # subcells is only inside this parent cell, not across to
        # other parent cells
        vtx_ids = {}
        for vtxs in subdiv_cells:
            subcell_id = len(subcells)

            # Add the vertices
            subcell = []
            for xN in vtxs:
                if xN not in vtx_ids:
                    vtx_ids[xN] = len(subvertices)
                    subvertices.append(xN)
                subcell.append(vtx_ids[xN])

            # Add the cell
            subcells.append(subcell)
            parent_cell[subcell_id] = cid

    # Create a new dolfin.Mesh
    subdivided_mesh = dolfin.Mesh()
    editor = dolfin.MeshEditor()
    editor.open(subdivided_mesh, 2, 2)

    # Add the vertices (nodes)
    editor.init_vertices(len(subvertices))
    for iv, vtx in enumerate(subvertices):
        editor.add_vertex(iv, vtx[0], vtx[1])

    # Add the cells (triangular elements)
    # Loop over the vertices and build two cells for each square
    # where the selected vertex is in the lower left corner
    editor.init_cells(len(subcells))
    for ic, subcell in enumerate(subcells):
        editor.add_cell(ic, subcell[0], subcell[1], subcell[2])

    # Close the mesh for editing
    editor.close()

    return subdivided_mesh, parent_cell


def subdivide_cell(n, x0, x1, x2):
    """
    Recursive refinement of triangular cells

    @param int n: number of refinmenents, n >= 0
    @param tuple(int) x0: vertex coordinate 0
    @param tuple(int) x1: vertex coordinate 1
    @param tuple(int) x2: vertex coordinate 2
    """
    if n == 0:
        return [(x0, x1, x2)]
    else:
        # Divide the cell such that the longest facet is halfed
        l0 = (x1[0] - x2[0]) ** 2 + (x1[1] - x2[1]) ** 2
        l1 = (x0[0] - x2[0]) ** 2 + (x0[1] - x2[1]) ** 2
        l2 = (x0[0] - x1[0]) ** 2 + (x0[1] - x1[1]) ** 2
        if l0 > max(l1, l2):
            xF = ((x1[0] + x2[0]) / 2, (x1[1] + x2[1]) / 2)
            return subdivide_cell(n - 1, x0, x1, xF) + subdivide_cell(n - 1, x0, xF, x2)
        elif l1 > max(l0, l2):
            xF = ((x0[0] + x2[0]) / 2, (x0[1] + x2[1]) / 2)
            return subdivide_cell(n - 1, x1, x2, xF) + subdivide_cell(n - 1, x1, xF, x0)
        else:
            xF = ((x0[0] + x1[0]) / 2, (x0[1] + x1[1]) / 2)
            return subdivide_cell(n - 1, x2, x0, xF) + subdivide_cell(n - 1, x2, xF, x1)
