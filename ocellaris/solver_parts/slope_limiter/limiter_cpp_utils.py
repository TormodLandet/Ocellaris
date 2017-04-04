import numpy
from ocellaris.cpp import load_module
cpp_mod = load_module('hierarchical_taylor')


class SlopeLimiterInput(object):
    def __init__(self, mesh, vertices, vertex_coordinates, num_neighbours,
                 neighbours, cell_dofs_V, cell_dofs_V0, limit_cell):
        """
        This class stores the connectivity and dof maps necessary to
        perform slope limiting in an efficient manner in the C++ code 
        """
        # Flatten 2D arrays for easy transfer to C++
        max_neighbours = neighbours.shape[1]
        flat_neighbours = neighbours.flatten()
        flat_cell_dofs = numpy.array(cell_dofs_V, dtype=numpy.intc).flatten()
        flat_cell_dofs_dg0 = numpy.array(cell_dofs_V0, dtype=numpy.intc).flatten()
        
        tdim = mesh.topology().dim()
        num_cells_owned = mesh.topology().ghost_offset(tdim)
        
        # Store coordinates for the three vertices plus the cell center for each cell
        stride = (3 + 1) * 2
        flat_vertex_coordinates = numpy.zeros(num_cells_owned * stride, float)
        for icell in xrange(num_cells_owned):
            cell_vertices = [vertex_coordinates[iv] for iv in vertices[icell]]
            center_pos_x = (cell_vertices[0][0] + cell_vertices[1][0] + cell_vertices[2][0]) / 3
            center_pos_y = (cell_vertices[0][1] + cell_vertices[1][1] + cell_vertices[2][1]) / 3
            istart = icell * stride
            flat_vertex_coordinates[istart+0:istart+2] = cell_vertices[0]
            flat_vertex_coordinates[istart+2:istart+4] = cell_vertices[1]
            flat_vertex_coordinates[istart+4:istart+6] = cell_vertices[2]
            flat_vertex_coordinates[istart+6:istart+8] = (center_pos_x, center_pos_y)

        # Call the C++ method that makes the arrays available to the C++ limiter
        self.cpp_obj = cpp_mod.SlopeLimiterInput()
        self.cpp_obj.set_arrays(num_cells_owned,
                                max_neighbours,
                                num_neighbours,
                                flat_neighbours,
                                flat_cell_dofs,
                                flat_cell_dofs_dg0,
                                flat_vertex_coordinates,
                                limit_cell)
    
    def set_global_bounds(self, global_min, global_max):
        """
        Set the minimum and maximum allowable field values for
        the limited field
        """
        self.global_min = global_min
        self.global_max = global_max
    
    def set_boundary_values(self, boundary_dof_type, boundary_dof_value):
        """
        Transfer the boundary dof data to the C++ side
        """
        assert boundary_dof_type.min() >= 0 and boundary_dof_type.max() <= 2
        self.cpp_obj.set_boundary_values(boundary_dof_type, boundary_dof_value)
