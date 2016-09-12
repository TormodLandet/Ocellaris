# enconding: utf8
from __future__ import division
import numpy
import dolfin as df
from ocellaris.cpp import load_module
from ocellaris.utils import ocellaris_error, verify_key
from . import register_slope_limiter, SlopeLimiterBase


@register_slope_limiter('nodal')
class SlopeLimiterNodal(SlopeLimiterBase):
    description = 'Ensures dof node values are not creating local extrema'
    
    def __init__(self, phi_name, phi, filter_method='nofilter', use_cpp=True):
        """
        Limit the slope of the given scalar to obtain boundedness
        """
        # Verify input
        V = phi.function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        loc = 'nodal slope limiter'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('function shape', phi.ufl_shape, [()], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        verify_key('filter', filter_method, ('nofilter', 'minmax'), loc)
        
        # Store input
        self.phi_name = phi_name
        self.phi = phi
        self.degree = V.ufl_element().degree()
        self.mesh = mesh
        self.filter = filter_method
        self.use_cpp = use_cpp
        
        # Exceedance is a secondary output of the limiter and is calculated
        # as the maximum correction performed for each cell 
        V0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.excedance = df.Function(V0)
        
        # Fast access to cell dofs
        dm, dm0 = V.dofmap(), V0.dofmap()
        indices = range(self.mesh.num_cells())
        cell_dofs_V = [tuple(dm.cell_dofs(i)) for i in indices]
        cell_dofs_V0 = [int(dm0.cell_dofs(i)) for i in indices]
        
        # Find the neighbour cells for each dof
        num_neighbours, neighbours = get_dof_neighbours(V)
        
        # Get the quadrature weights for each dof in each cell
        # (needed to compute the average value)
        quadrature_weights = get_quadrature_weights(V)
        
        # Get indices for get_local on the DGX vector
        im = dm.index_map()
        n_local_and_ghosts = im.size(im.MapSize_ALL)
        intc = numpy.intc
        self.local_indices_dgX = numpy.arange(n_local_and_ghosts, dtype=intc)
        
        # Flatten 2D arrays for easy transfer to C++
        self.num_neighbours = num_neighbours
        self.max_neighbours = neighbours.shape[1]
        self.num_cell_dofs = quadrature_weights.shape[1]
        self.flat_neighbours = neighbours.flatten()
        self.flat_cell_dofs = numpy.array(cell_dofs_V, dtype=intc).flatten()
        self.flat_cell_dofs_dg0 = numpy.array(cell_dofs_V0, dtype=intc).flatten()
        self.flat_weights = quadrature_weights.flatten()
        self.cpp_mod = load_module('slope_limiter_basic')
        
        # The initial maximum and minimum values in the boundeness filter are cached
        self._filter_cache = None
    
    def run(self):
        """
        Perform the slope limiting and filtering
        """
        assert self.degree == 1, 'Slope limiter only implemented for linear elements'
        
        # Get local values + the ghost values
        results = numpy.zeros(len(self.local_indices_dgX), float)
        self.phi.vector().get_local(results, self.local_indices_dgX)
        
        # Get the limiter implementation based on polynomial degree
        # and the implementation language (Python prototype or the
        # default faster C++ implementation)
        if self.use_cpp:
            if self.degree == 1:
                limiter = self.cpp_mod.slope_limiter_basic_dg1
            else:
                ocellaris_error('Slope limiter error',
                                'C++ slope limiter does not support degree %d' % self.degree)
        else:
            if self.degree == 1:
                limiter = slope_limiter_basic_dg1
            else:
                ocellaris_error('Slope limiter error',
                                'Python slope limiter does not support degree %d' % self.degree)
        
        num_cells_all = self.mesh.num_cells()
        tdim = self.mesh.topology().dim()
        num_cells_owned = self.mesh.topology().ghost_offset(tdim)
        exceedances = self.excedance.vector().get_local()
        
        limiter(self.num_neighbours,
                num_cells_all,
                num_cells_owned,
                self.num_cell_dofs,
                self.max_neighbours,
                self.flat_neighbours,
                self.flat_cell_dofs,
                self.flat_cell_dofs_dg0,
                self.flat_weights,
                exceedances,
                results)
        
        self.excedance.vector().set_local(exceedances)
        self.excedance.vector().apply('insert')
        
        # Run post processing filter
        if self.filter == 'minmax':
            self._run_minmax_filter(results)
        
        self.phi.vector().set_local(results, self.local_indices_dgX)
        self.phi.vector().apply('insert')
    
    def _run_minmax_filter(self, results):
        """
        Make sure the DOF values are inside the min and max values from the
        first time step (enforce boundedness)
        """
        if self._filter_cache is None:
            mi = df.MPI.min(df.mpi_comm_world(), float(results.min()))
            ma = df.MPI.max(df.mpi_comm_world(), float(results.max()))
            self._filter_cache = mi, ma
            return
        minval, maxval = self._filter_cache
        numpy.clip(results, minval, maxval, results)


def get_dof_neighbours(V):
    """
    Given a DG function space find, for each dof, the indices
    of the cells with dofs at the same locations
    """
    dm = V.dofmap()
    gdim = V.mesh().geometry().dim()
    num_cells_all = V.mesh().num_cells()
    dof_coordinates = V.tabulate_dof_coordinates().reshape((-1, gdim))
    
    # Get "owning cell" indices for all dofs
    cell_for_dof = [None] * V.dim()
    for ic in xrange(num_cells_all):
        dofs = dm.cell_dofs(ic)
        for dof in dofs:
            assert cell_for_dof[dof] is None
            cell_for_dof[dof] = ic
    
    # Map dof coordinate to dofs, this is for DG so multiple dofs
    # will share the same location
    coord_to_dofs = {}
    max_neighbours = 0
    for dof in xrange(len(dof_coordinates)):
        coord = tuple(round(x, 5) for x in dof_coordinates[dof])
        dofs = coord_to_dofs.setdefault(coord, [])
        dofs.append(dof)
        max_neighbours = max(max_neighbours, len(dofs)-1)
    
    # Find number of neighbour cells and their indices for each dof
    num_neighbours = numpy.zeros(V.dim(), numpy.intc)
    neighbours = numpy.zeros((V.dim(), max_neighbours), numpy.intc) - 1
    for nbs in coord_to_dofs.values():
        # Loop through dofs at this location
        for dof in nbs:
            # Loop through the dofs neighbours
            for nb in nbs:
                # Skip the dof itself
                if dof == nb:
                    continue
                # Get the neighbours "owning cell" index and store this
                nb_cell = cell_for_dof[nb]
                nn_prev = num_neighbours[dof]
                neighbours[dof,nn_prev] = nb_cell
                num_neighbours[dof] += 1
    
    return num_neighbours, neighbours


def get_quadrature_weights(V):
    """
    Get quadrature weights needed to compute the cell average
    of the functions given the dof values
    """
    degree = V.ufl_element().degree()
    N = V.mesh().num_cells()
    
    if degree == 1:
        weights = numpy.ones((N, 3), float)/3 
    else:
        ocellaris_error('Quadrature weight error',
                        'Cannot compute quadrature weights on degree %d function' % degree)
        
    return weights


###################################################################################################
# Python implementations of nodal slope limiters
#
# The interface is a bit un-Pythonic, but this is done to be able to have the
# same interface for the C++ and the Python implementations. The Python
# implementations are meant to be prototypes and QA of the C++ implementations


def slope_limiter_basic_dg1(num_neighbours, num_cells_all, num_cells_owned, num_cell_dofs,
                            max_neighbours, flat_neighbours, flat_cell_dofs, flat_cell_dofs_dg0,
                            flat_weights, exceedances, results):
    """
    Perform nodal slope limiting of a DG1 function
    """
    # Cell averages
    averages = numpy.zeros(num_cells_all, float)
    for ic in range(num_cells_all):
        dofs = flat_cell_dofs[ic * 3: (ic + 1)*3]
        vals = [results[dof] for dof in dofs]
        weights = flat_weights[ic * 3: (ic + 1)*3]
        averages[ic] = numpy.dot(vals, weights)
    
    # Modify dof values
    dof_range = range(3)
    for ic in range(num_cells_owned):
        dofs = flat_cell_dofs[ic * 3: (ic + 1)*3]
        vals = [results[dof] for dof in dofs]
        weights = flat_weights[ic * 3: (ic + 1)*3]
        avg = sum(vals) / 3
        
        excedance = 0
        for idof in dof_range:
            dof = dofs[idof]
            n_nbs = num_neighbours[dof]
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
        
        exceedances[flat_cell_dofs_dg0[ic]] = excedance
        if excedance == 0:
            continue
        
        # Modify the results to limit the slope
        
        # Find the new average and which vertices can be adjusted to obtain the correct average
        new_avg = numpy.dot(vals, weights)
        eps = 0
        moddable = [0]*len(dof_range)
        if abs(avg - new_avg) > 1e-15:
            if new_avg > avg:
                for idof in dof_range:
                    if vals[idof] > avg:
                        moddable[idof] = 1
            else:
                for idof in dof_range:
                    if vals[idof] < avg:
                        moddable[idof] = 1
            
            # Get number of vertex values that can be modified and catch
            # possible floating point problems with the above comparisons
            nmod = sum(moddable)
            
            if nmod == 0:
                assert abs(excedance) < 1e-14, 'Nmod=0 with exceedance=%r' % excedance
            else:
                eps = (avg - new_avg) * 3 / nmod
        
        # Modify the vertices to obtain the correct average
        for idof in dof_range:
            results[dofs[idof]] = vals[idof] + eps * moddable[idof]
    
    return exceedances
