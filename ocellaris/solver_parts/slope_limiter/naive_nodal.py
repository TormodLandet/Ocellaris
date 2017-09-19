import numpy
import dolfin as df
from ocellaris.cpp import load_module
from ocellaris.utils import ocellaris_error, verify_key, get_dof_neighbours
from . import register_slope_limiter, SlopeLimiterBase


@register_slope_limiter('NaiveNodal')
class NaiveNodalSlopeLimiter(SlopeLimiterBase):
    description = 'Ensures dof node values are not themselves a local extrema'
    
    def __init__(self, phi_name, phi, boundary_condition, output_name=None, use_cpp=True, enforce_bounds=False):
        """
        Limit the slope of the given scalar to obtain boundedness
        """
        # Verify input
        V = phi.function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        loc = 'NaiveNodal slope limiter'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (0, 1, 2), loc)
        verify_key('function shape', phi.ufl_shape, [()], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        
        # Store input
        self.phi_name = phi_name
        self.phi = phi
        self.degree = degree
        self.mesh = mesh
        self.use_cpp = use_cpp
        self.enforce_global_bounds = enforce_bounds
        
        if output_name is None:
            output_name = phi_name
        
        # Exceedance is a secondary output of the limiter and is calculated
        # as the maximum correction performed for each cell 
        V0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.exceedance = df.Function(V0)
        name = 'SlopeLimiterExceedance_%s' % output_name
        self.exceedance.rename(name, name)
        self.additional_plot_funcs = [self.exceedance]
        
        # No limiter needed for piecewice constant functions
        if degree == 0:
            return
        
        # Fast access to cell dofs
        dm, dm0 = V.dofmap(), V0.dofmap()
        indices = list(range(self.mesh.num_cells()))
        cell_dofs_V = [tuple(dm.cell_dofs(i)) for i in indices]
        cell_dofs_V0 = [int(dm0.cell_dofs(i)) for i in indices]
        
        # Find the neighbour cells for each dof
        num_neighbours, neighbours = get_dof_neighbours(V)
        
        # Remove boundary dofs from limiter
        num_neighbours[boundary_condition != 0] = 0
        
        # Get indices for get_local on the DGX vector
        im = dm.index_map()
        n_local_and_ghosts = im.size(im.MapSize_ALL)
        intc = numpy.intc
        self.local_indices_dgX = numpy.arange(n_local_and_ghosts, dtype=intc)
        
        # Flatten 2D arrays for easy transfer to C++
        self.num_neighbours = num_neighbours
        self.max_neighbours = neighbours.shape[1]
        self.num_cell_dofs = 3 if self.degree == 1 else 6
        self.flat_neighbours = neighbours.flatten()
        self.flat_cell_dofs = numpy.array(cell_dofs_V, dtype=intc).flatten()
        self.flat_cell_dofs_dg0 = numpy.array(cell_dofs_V0, dtype=intc).flatten()
        self.cpp_mod = load_module('naive_nodal')
    
    def run(self):
        """
        Perform the slope limiting
        """
        # No limiter needed for piecewice constant functions
        if self.degree == 0:
            return
        
        # Get local values + the ghost values
        results = numpy.zeros(len(self.local_indices_dgX), float)
        self.phi.vector().get_local(results, self.local_indices_dgX)
        
        # Get the limiter implementation based on polynomial degree
        # and the implementation language (Python prototype or the
        # default faster C++ implementation)
        if self.use_cpp:
            if self.degree in (1,):
                limiter = self.cpp_mod.naive_nodal_slope_limiter_dg1
            else:
                ocellaris_error('NaiveNodal slope limiter error',
                                'C++ slope limiter does not support degree %d' % self.degree)
        else:
            if self.degree in (1, 2):
                limiter = slope_limiter_nodal_dg
            else:
                ocellaris_error('NaiveNodal slope limiter error',
                                'Python slope limiter does not support degree %d' % self.degree)
        
        num_cells_all = self.mesh.num_cells()
        tdim = self.mesh.topology().dim()
        num_cells_owned = self.mesh.topology().ghost_offset(tdim)
        exceedances = self.exceedance.vector().get_local()
        
        limiter(self.num_neighbours,
                num_cells_all,
                num_cells_owned,
                self.num_cell_dofs,
                self.max_neighbours,
                self.flat_neighbours,
                self.flat_cell_dofs,
                self.flat_cell_dofs_dg0,
                exceedances,
                results)
        
        self.exceedance.vector().set_local(exceedances)
        self.exceedance.vector().apply('insert')
        
        # Enforce boundedness by cliping values to the bounds
        if self.enforce_global_bounds:
            minval, maxval = self.global_bounds
            numpy.clip(results, minval, maxval, results)
        
        self.phi.vector().set_local(results, self.local_indices_dgX)
        self.phi.vector().apply('insert')    


###################################################################################################
# Python implementations of nodal slope limiters
#
# The interface is a bit un-Pythonic, but this is done to be able to have the
# same interface for the C++ and the Python implementations. The Python
# implementations are meant to be prototypes and QA of the C++ implementations


def slope_limiter_nodal_dg(num_neighbours, num_cells_all, num_cells_owned, num_cell_dofs,
                           max_neighbours, flat_neighbours, flat_cell_dofs, flat_cell_dofs_dg0,
                           exceedances, results):
    """
    Perform nodal slope limiting of a DG1 function. Also works for DG2
    functions, but since the local minimum or maximum may be in between
    dof locations the slope limits are not as strict as they should be
    since this implementation only limits dof values. 
    """
    C = num_cell_dofs
    
    # Compute cell averages before any modification
    averages = numpy.zeros(num_cells_all, float)
    for ic in range(num_cells_all):
        dofs = flat_cell_dofs[ic * C: (ic + 1)*C]
        vals = [results[dof] for dof in dofs]
        averages[ic] = sum(vals[-3:]) / 3
    
    # Modify dof values
    dof_range = list(range(C))
    dof_range_modable = dof_range[-3:]
    for ic in range(num_cells_owned):
        dofs = flat_cell_dofs[ic * C: (ic + 1)*C]
        vals = [results[dof] for dof in dofs]
        avg = averages[ic]
        
        excedance = 0
        cell_on_boundary = False
        for idof in dof_range:
            dof = dofs[idof]
            n_nbs = num_neighbours[dof]
            
            if n_nbs == 0:
                # Skip boundary dofs which are marked with zero neighbours
                cell_on_boundary = True
                break
            
            nbs = flat_neighbours[dof * max_neighbours: dof * max_neighbours + n_nbs]
            nb_vals = [averages[nb] for nb in nbs]
            
            # Find highest and lowest value in the connected cells
            lo = hi = avg
            for cell_avg in nb_vals:
                lo = min(lo, cell_avg)
                hi = max(hi, cell_avg)
            
            vtx_val = vals[idof]
            if vtx_val < lo:
                vals[idof] = lo
                ex = vtx_val - lo
                if abs(excedance) < abs(ex):
                    excedance = ex
            elif vtx_val > hi:
                vals[idof] = hi
                ex = vtx_val - hi
                if abs(excedance) < abs(ex):
                    excedance = ex
        
        if cell_on_boundary:
            continue
        
        exceedances[flat_cell_dofs_dg0[ic]] = excedance
        if excedance == 0:
            continue
        
        # Modify the results to limit the slope
        
        # Find the new average and which vertices can be adjusted to obtain the correct average
        new_avg = sum(vals[-3:]) / 3
        eps = 0
        moddable = [0]*len(dof_range)
        if abs(avg - new_avg) > 1e-15:
            if new_avg > avg:
                for idof in dof_range_modable:
                    if vals[idof] > avg:
                        moddable[idof] = 1
            else:
                for idof in dof_range_modable:
                    if vals[idof] < avg:
                        moddable[idof] = 1
            
            # Get number of vertex values that can be modified and catch
            # possible floating point problems with the above comparisons
            nmod = sum(moddable)
            assert nmod <= 3
            
            if nmod == 0:
                assert abs(excedance) < 1e-14, 'Nmod=0 with exceedance=%r' % excedance
            else:
                eps = (avg - new_avg) * 3 / nmod 
        
        # Modify the vertices to obtain the correct average
        for idof, dof in enumerate(dofs):
            results[dof] = vals[idof] + eps * moddable[idof]
        
        # Verify that the average was not changed
        vals_final = [results[dof] for dof in dofs]
        final_avg = sum(vals_final[-3:]) / 3
        error = abs(averages[ic] - final_avg)/averages[ic]
        if error > 1e-10:
            new_vals = [vals[i] + eps*moddable[i] for i in dof_range]
            print(num_cell_dofs, repr(error))
            print(repr(avg), repr(new_avg))
            print(repr(avg - new_avg), repr(eps), 'sssssssssssssssssssssssssssss')
            print(dofs, dof_range, dof_range_modable)
            print('%r %r %r' % (averages[ic], new_avg, final_avg))
            print(moddable, eps, nmod)
            print(repr(averages[ic] - (sum(vals[-3:]) + eps*nmod)/3))
            print(vals)
            print(new_vals)
            print([results[dof] for dof in dofs])
            print(repr(averages[ic] - sum(new_vals) / 3))
            print('--------------------------------------------------------------------------')
            exit()
        assert error < 1e-12, 'Got large difference in old and new average: %r' % error 
    
    return exceedances
