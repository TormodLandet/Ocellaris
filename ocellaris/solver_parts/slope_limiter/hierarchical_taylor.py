# encoding: utf8
from __future__ import division
import numpy
import dolfin as df
from ocellaris.cpp import load_module
from ocellaris.utils import verify_key, get_dof_neighbours
from ocellaris.utils import lagrange_to_taylor, taylor_to_lagrange
from . import register_slope_limiter, SlopeLimiterBase


@register_slope_limiter('HierarchicalTaylor')
class HierarchicalTaylorSlopeLimiter(SlopeLimiterBase):
    description = 'Uses a Taylor DG decomposition to limit derivatives at the vertices in a hierarchical manner'
    
    def __init__(self, phi_name, phi, boundary_condition=None, filter_method='nofilter',
                 use_cpp=True, output_name=None):
        """
        Limit the slope of the given scalar to obtain boundedness
        """
        # Verify input
        V = phi.function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        loc = 'HierarchalTaylor slope limiter'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (0, 1, 2), loc)
        verify_key('function shape', phi.ufl_shape, [()], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        verify_key('filter', filter_method, ('nofilter',), loc)
        
        # Store input
        self.phi_name = phi_name
        self.phi = phi
        self.degree = degree
        self.mesh = mesh
        self.filter = filter_method
        self.use_cpp = use_cpp
        
        if output_name is None:
            output_name = phi_name
        
        # Alpha factors are secondary outputs
        V0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.alpha_funcs = []
        for i in range(degree):
            func = df.Function(V0)
            name = 'SlopeLimiterAlpha%d_%s' % (i+1, output_name)
            func.rename(name, name)
            self.alpha_funcs.append(func)
        self.additional_plot_funcs = self.alpha_funcs
        
        # Intermediate DG Taylor function space
        self.taylor = df.Function(V)
        
        # No limiter needed for piecewice constant functions
        if degree == 0:
            return
        
        # Find the neighbour cells for each dof
        num_neighbours, neighbours = get_dof_neighbours(V)
        self.num_neighbours = num_neighbours
        self.neighbours = neighbours
        
        # Remove boundary dofs from limiter
        if boundary_condition is not None:
            num_neighbours[boundary_condition != 0] = 0
        
        # Fast access to cell dofs
        dm, dm0 = V.dofmap(), V0.dofmap()
        indices = range(self.mesh.num_cells())
        self.cell_dofs_V = [tuple(dm.cell_dofs(i)) for i in indices]
        self.cell_dofs_V0 = [int(dm0.cell_dofs(i)) for i in indices]
        
        # Find vertices for each cell
        mesh.init(2, 0)
        connectivity_CV = mesh.topology()(2, 0)
        vertices = []
        for ic in range(self.mesh.num_cells()):
            vnbs = tuple(connectivity_CV(ic))
            vertices.append(vnbs)
        self.vertices = vertices
        self.vertex_coordinates = mesh.coordinates()
        
        if use_cpp:
            # Flatten 2D arrays for easy transfer to C++
            self.num_neighbours = num_neighbours
            self.max_neighbours = neighbours.shape[1]
            self.flat_neighbours = neighbours.flatten()
            self.flat_cell_dofs = numpy.array(self.cell_dofs_V, dtype=numpy.intc).flatten()
            self.flat_cell_dofs_dg0 = numpy.array(self.cell_dofs_V0, dtype=numpy.intc).flatten()
            self.cpp_mod = load_module('hierarchical_taylor')
            
            tdim = self.mesh.topology().dim()
            self.num_cells_owned = self.mesh.topology().ghost_offset(tdim)
            
            # Store coordinates for the three vertices plus the cell center for each cell
            stride = (3 + 1) * 2
            self.flat_vertex_coordinates = numpy.zeros(self.num_cells_owned * stride, float)
            for icell in xrange(self.num_cells_owned):
                cell_vertices = [self.vertex_coordinates[iv] for iv in self.vertices[icell]]
                center_pos_x = (cell_vertices[0][0] + cell_vertices[1][0] + cell_vertices[2][0]) / 3
                center_pos_y = (cell_vertices[0][1] + cell_vertices[1][1] + cell_vertices[2][1]) / 3
                istart = icell * stride
                self.flat_vertex_coordinates[istart+0:istart+2] = cell_vertices[0]
                self.flat_vertex_coordinates[istart+2:istart+4] = cell_vertices[1]
                self.flat_vertex_coordinates[istart+4:istart+6] = cell_vertices[2]
                self.flat_vertex_coordinates[istart+6:istart+8] = (center_pos_x, center_pos_y)
    
    def run(self):
        """
        Perform slope limiting of DG Lagrange functions
        """
        timer = df.Timer('Ocellaris HierarchalTaylorSlopeLimiter')
        
        # Update the Taylor function with the current Lagrange values
        lagrange_to_taylor(self.phi, self.taylor)
        taylor_arr = self.taylor.vector().get_local()
        alpha_arrs = [alpha.vector().get_local() for alpha in self.alpha_funcs]
        
        if self.use_cpp:
            self._run_cpp(taylor_arr, alpha_arrs)
        elif self.degree == 1:
            self._run_dg1(taylor_arr, alpha_arrs[0])
        elif self.degree == 2:
            self._run_dg2(taylor_arr, alpha_arrs[0], alpha_arrs[1])
            
        # Update the Lagrange function with the limited Taylor values
        self.taylor.vector().set_local(taylor_arr)
        self.taylor.vector().apply('insert')
        taylor_to_lagrange(self.taylor, self.phi)
        
        # Update the secondary output arrays, alphas
        for alpha, alpha_arr in zip(self.alpha_funcs, alpha_arrs):
            alpha.vector().set_local(alpha_arr)
            alpha.vector().apply('insert')
        
        timer.stop()
        
    def _run_cpp(self, taylor_arr, alpha_arrs):
        if self.degree == 1:
            limiter = self.cpp_mod.hierarchical_taylor_slope_limiter_dg1
        elif self.degree == 2:
            limiter = self.cpp_mod.hierarchical_taylor_slope_limiter_dg2
        
        limiter(self.num_neighbours,
                self.num_cells_owned,
                self.max_neighbours,
                self.flat_neighbours,
                self.flat_cell_dofs,
                self.flat_cell_dofs_dg0,
                self.flat_vertex_coordinates,
                taylor_arr,
                *alpha_arrs)
    
    def _run_dg1(self, taylor_arr, alpha_arr):
        """
        Perform slope limiting of a DG1 function
        """
        lagrange_arr = self.phi.vector().get_local()
        
        V = self.phi.function_space()
        mesh = V.mesh()
        tdim = mesh.topology().dim()
        num_cells_owned = mesh.topology().ghost_offset(tdim)
        
        for icell in xrange(num_cells_owned):
            dofs = self.cell_dofs_V[icell]
            center_value = taylor_arr[dofs[0]]
            
            # Find the minimum slope limiter coefficient alpha 
            alpha = 1.0
            for i in xrange(3):
                dof = dofs[i]
                
                # Find vertex neighbours minimum and maximum values
                minval = maxval = center_value
                for nb in self.neighbours[dof]:
                    nb_center_val_dof = self.cell_dofs_V[nb][0]
                    nb_val = taylor_arr[nb_center_val_dof]
                    minval = min(minval, nb_val)
                    maxval = max(maxval, nb_val)
                
                vertex_value = lagrange_arr[dof]
                if vertex_value > center_value:
                    alpha = min(alpha, (maxval - center_value)/(vertex_value - center_value))
                elif vertex_value < center_value:
                    alpha = min(alpha, (minval - center_value)/(vertex_value - center_value))
            
            alpha_arr[self.cell_dofs_V0[icell]] = alpha
            taylor_arr[dofs[1]] *= alpha
            taylor_arr[dofs[2]] *= alpha
    
    def _run_dg2(self, taylor_arr, alpha1_arr, alpha2_arr):
        """
        Perform slope limiting of a DG2 function
        """
        V = self.phi.function_space()
        mesh = V.mesh()
        tdim = mesh.topology().dim()
        num_cells_owned = mesh.topology().ghost_offset(tdim)
        
        # Slope limit one cell at a time
        for icell in xrange(num_cells_owned):
            dofs = self.cell_dofs_V[icell]
            assert len(dofs) == 6
            center_values = [taylor_arr[dof] for dof in dofs]
            (center_phi, center_phix, center_phiy, center_phixx, 
                center_phiyy, center_phixy) = center_values
            
            cell_vertices = [self.vertex_coordinates[iv] for iv in self.vertices[icell]]
            center_pos_x = (cell_vertices[0][0] + cell_vertices[1][0] + cell_vertices[2][0]) / 3
            center_pos_y = (cell_vertices[0][1] + cell_vertices[1][1] + cell_vertices[2][1]) / 3
            assert len(cell_vertices) == 3
            
            # Find the minimum slope limiter coefficient alpha of the φ, dφdx and dφ/dy terms
            alpha = [1.0] * 3
            for taylor_dof in (0, 1, 2): 
                for ivert in xrange(3):
                    dof = dofs[ivert]
                    dx = cell_vertices[ivert][0] - center_pos_x
                    dy = cell_vertices[ivert][1] - center_pos_y
                    
                    # Find vertex neighbours minimum and maximum values
                    base_value = center_values[taylor_dof]
                    minval = maxval = base_value
                    for nb in self.neighbours[dof]:
                        nb_center_val_dof = self.cell_dofs_V[nb][taylor_dof]
                        nb_val = taylor_arr[nb_center_val_dof]
                        minval = min(minval, nb_val)
                        maxval = max(maxval, nb_val)
                    
                    # Compute vertex value
                    if taylor_dof == 0:
                        # Function value at the vertex (linear reconstruction)
                        vertex_value = center_phi + center_phix * dx + center_phiy * dy
                    elif taylor_dof == 1:
                        # Derivative in x direction at the vertex  (linear reconstruction)
                        vertex_value = center_phix + center_phixx * dx + center_phixy * dy
                    else:
                        # Derivative in y direction at the vertex  (linear reconstruction)
                        vertex_value = center_phiy + center_phiyy * dy + center_phixy * dx
                    
                    # Compute the slope limiter coefficient alpha
                    if vertex_value > base_value:
                        a = (maxval - base_value) / (vertex_value - base_value)
                    elif vertex_value < base_value:
                        a = (minval - base_value) / (vertex_value - base_value)
                    else:
                        a = 1
                    alpha[taylor_dof] = min(alpha[taylor_dof], a)
            
            alpha2 = min(alpha[1], alpha[2])
            alpha1 = max(alpha[0], alpha2)
            
            taylor_arr[dofs[1]] *= alpha1
            taylor_arr[dofs[2]] *= alpha1
            taylor_arr[dofs[3]] *= alpha2
            taylor_arr[dofs[4]] *= alpha2
            taylor_arr[dofs[5]] *= alpha2
            
            dof_dg0 = self.cell_dofs_V0[icell] 
            alpha1_arr[dof_dg0] = alpha1
            alpha2_arr[dof_dg0] = alpha2
