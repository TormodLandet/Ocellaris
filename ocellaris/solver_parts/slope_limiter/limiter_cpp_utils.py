# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
from ocellaris.cpp import load_module
from ocellaris.utils import get_dof_neighbours


class SlopeLimiterInput(object):
    def __init__(self, mesh, V, V0, use_cpp=True, trust_robin_dval=True):
        """
        This class stores the connectivity and dof maps necessary to
        perform slope limiting in an efficient manner in the C++ code
        """
        # Find the neighbour cells for each dof
        num_neighbours, neighbours = get_dof_neighbours(V)
        self.num_neighbours = num_neighbours
        self.neighbours = neighbours

        # Fast access to cell dofs
        dm, dm0 = V.dofmap(), V0.dofmap()
        indices = list(range(mesh.num_cells()))
        self.cell_dofs_V = numpy.array([dm.cell_dofs(i) for i in indices], dtype=numpy.intc)
        self.cell_dofs_V0 = numpy.array([dm0.cell_dofs(i) for i in indices], dtype=numpy.intc)

        tdim = mesh.topology().dim()
        self.num_cells_owned = mesh.topology().ghost_offset(tdim)

        # Find vertices for each cell
        mesh.init(tdim, 0)
        connectivity_CV = mesh.topology()(tdim, 0)
        cell_vertices = []
        for ic in range(mesh.num_cells()):
            cell_vertices.append(connectivity_CV(ic))
        self.cell_vertices = numpy.array(cell_vertices, dtype=numpy.intc)
        self.vertex_coordinates = mesh.coordinates()
        assert self.vertex_coordinates.shape[1] == tdim

        # Call the C++ method that makes the arrays available to the C++ limiter
        if use_cpp:
            self.cpp_mod = load_module('hierarchical_taylor')
            self.cpp_obj = self.cpp_mod.SlopeLimiterInput()
            self.cpp_obj.set_arrays(
                self.num_cells_owned,
                num_neighbours,
                neighbours,
                self.cell_dofs_V,
                self.cell_dofs_V0,
                self.cell_vertices,
                self.vertex_coordinates,
            )
            self.cpp_obj.trust_robin_dval = trust_robin_dval

    def set_global_bounds(self, global_min, global_max):
        """
        Set the minimum and maximum allowable field values for
        the limited field
        """
        self.cpp_obj.global_min = global_min
        self.cpp_obj.global_max = global_max

    def set_limit_cell(self, limit_cell):
        """
        Decide whether each cell should be limited. Accepts an array of length
        num_cells_owned which is mesh.topology().ghost_offset(mesh.topology.dim())

        There reason for int instead of bool is purely ease of interfacing with
        the C++ SWIG wrapper

        @param numpy.ndarray(int) limit_cell: 1 for cells that should be limited, else 0
        """
        self.cpp_obj.set_limit_cell(limit_cell)

    def set_boundary_values(self, boundary_dof_type, boundary_dof_value, enforce):
        """
        Transfer the boundary dof data to the C++ side
        """
        assert boundary_dof_type.min() >= 0 and boundary_dof_type.max() <= 4
        self.cpp_obj.set_boundary_values(boundary_dof_type, boundary_dof_value, enforce)

    def get_cpp_mod(self):
        return self.cpp_mod
