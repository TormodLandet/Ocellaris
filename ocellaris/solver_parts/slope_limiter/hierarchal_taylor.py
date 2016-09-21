# enconding: utf8
from __future__ import division
import numpy
import dolfin as df
from ocellaris.cpp import load_module
from ocellaris.utils import ocellaris_error, verify_key, get_dof_neighbours
from ocellaris.utils import lagrange_to_taylor, taylor_to_lagrange
from . import register_slope_limiter, SlopeLimiterBase


@register_slope_limiter('HierarchalTaylor')
class HierarchalTaylorSlopeLimiter(SlopeLimiterBase):
    description = 'Uses a Taylor DG decomposition to limit derivatives at the vertices in a hierarchal manner'
    
    def __init__(self, phi_name, phi, boundary_condition, filter_method='nofilter', use_cpp=True):
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
        verify_key('slope limited degree', degree, (0, 1), loc)
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
        
        # Exceedance is a secondary output of the limiter and is calculated
        # as the minimum slope limiter coefficient alpha for each cell
        V0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.excedance = df.Function(V0)
        
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
    
    def run(self):
        """
        Perform slope limiting
        """
        # No limiter needed for piecewice constant functions
        if self.degree == 0:
            return
        
        # Update the Taylor function space with the new DG values
        lagrange_to_taylor(self.phi, self.taylor)
        lagrange_vals = self.phi.vector().get_local()
        taylor_vals = self.taylor.vector().get_local()
        alphas = self.excedance.vector().get_local()
        
        V = self.phi.function_space()
        mesh = V.mesh()
        tdim = mesh.topology().dim()
        num_cells_owned = mesh.topology().ghost_offset(tdim)
        
        for icell in xrange(num_cells_owned):
            dofs = self.cell_dofs_V[icell]
            center_value = taylor_vals[dofs[0]]
            
            # Find the minimum slope limiter coefficient alpha 
            alpha = 1.0
            for i in xrange(3):
                dof = dofs[i]
                if not self.num_neighbours[dof]:
                    alpha = 1.0
                    break
                
                # Find vertex neighbours minimum and maximum values
                minval = maxval = center_value
                for nb in self.neighbours[dof]:
                    nb_center_val_dof = self.cell_dofs_V[nb][0]
                    nb_val = taylor_vals[nb_center_val_dof]
                    minval = min(minval, nb_val)
                    maxval = max(maxval, nb_val)
                
                vertex_value = lagrange_vals[dof]
                if vertex_value > center_value:
                    alpha = min(alpha, (maxval - center_value)/(vertex_value - center_value))
                elif vertex_value < center_value:
                    alpha = min(alpha, (minval - center_value)/(vertex_value - center_value))
            
            alphas[self.cell_dofs_V0[icell]] = alpha
            taylor_vals[dofs[1]] *= alpha
            taylor_vals[dofs[2]] *= alpha
        
        # Update the DG Lagrange function space with the limited DG Taylor values
        self.taylor.vector().set_local(taylor_vals)
        taylor_to_lagrange(self.taylor, self.phi)
        
        self.excedance.vector().set_local(alphas)
