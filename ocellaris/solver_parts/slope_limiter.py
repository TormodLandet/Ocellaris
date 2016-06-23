# enconding: utf8
from __future__ import division
import numpy
import dolfin as df
from dolfin import dot, dx
from ocellaris.cpp import load_module
from ocellaris.utils import verify_key


LIMITER = 'none'
FILTER = 'nofilter'
USE_CPP = True


class SlopeLimiterBasic(object):
    def __init__(self, phi, filter_method='nofilter', use_cpp=True):
        """
        Limit the slope of the given scalar to obtain boundedness
        """
        self.phi = phi
        V = phi.function_space()
        self.degree = V.ufl_element().degree()
        self.mesh = V.mesh()
        
        assert V.ufl_element().family() == 'Discontinuous Lagrange'
        assert filter_method in ('nofilter', 'minmax')
        assert phi.ufl_shape == ()
        assert self.mesh.topology().dim() == 2
        
        self.filter = filter_method
        
        # Intermediate function spaces used in the limiter
        self.DG0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.DG1 = df.FunctionSpace(self.mesh, 'DG', 1)
        self.excedance = df.Function(self.DG0)
        
        # Projection to DG1 with pre factorized matrices for quick application
        if self.degree > 1:
            u, v = df.TrialFunction(self.DG1), df.TestFunction(self.DG1)
            a = dot(u, v)*dx
            L = dot(phi, v)*dx
            self.local_solver1 = df.LocalSolver(a, L)
            self.local_solver1.factorize()
        
        # Fast access to cell dofs
        dm0, dm1 = self.DG0.dofmap(), self.DG1.dofmap()
        indices = range(self.mesh.num_cells())
        self.cell_dofs_DG0 = list(int(dm0.cell_dofs(i)) for i in indices)
        self.cell_dofs_DG1 = list(tuple(dm1.cell_dofs(i)) for i in indices)
        
        # Connectivity from cell to vertex and vice versa
        self.mesh.init(2, 0)
        self.mesh.init(0, 2)
        connectivity_CV = self.mesh.topology()(2, 0)
        connectivity_VC = self.mesh.topology()(0, 2)
        
        # Find the neighbours cells for each vertex
        self.neighbours = []
        for iv in range(self.mesh.num_vertices()):
            cnbs = tuple(connectivity_VC(iv))
            self.neighbours.append(cnbs)
        
        # Find vertices for each cell
        self.vertices = []
        for ic in range(self.mesh.num_cells()):
            vnbs = tuple(connectivity_CV(ic))
            self.vertices.append(vnbs)
        
        self._filter_cache = None
        
        # Get indices for get_local on the DG1 vector
        im = dm1.index_map()
        n_local_and_ghosts = im.size(im.MapSize_ALL)
        intc = numpy.intc
        self.local_indices_dg1 = numpy.arange(n_local_and_ghosts, dtype=intc)
        
        self.use_cpp = use_cpp
        if use_cpp:
            self.num_neighbours = numpy.array([len(nbs) for nbs in self.neighbours], dtype=intc)
            self.max_neighbours = Nnmax = self.num_neighbours.max()
            self.flat_neighbours = numpy.zeros(len(self.neighbours)*Nnmax, dtype=intc) - 1
            for i, nbs in enumerate(self.neighbours):
                self.flat_neighbours[i*Nnmax:i*Nnmax + self.num_neighbours[i]] = self.neighbours[i]
            self.flat_cell_dofs = numpy.array(self.cell_dofs_DG1, dtype=intc).flatten()
            self.flat_cell_dofs_dg0 = numpy.array(self.cell_dofs_DG0, dtype=intc).flatten()
            self.flat_vertices = numpy.array(self.vertices, dtype=intc).flatten()
            self.cpp_mod = load_module('slope_limiter_basic')
    
    def run(self):
        """
        Perform the slope limiting and filtering
        """
        assert self.degree == 1, 'Slope limiter only implemented for linear elements'
        
        # Get local values + the ghost values
        results = numpy.zeros(len(self.local_indices_dg1), float)
        self.phi.vector().get_local(results, self.local_indices_dg1)
        
        if self.degree > 0:
            if self.use_cpp:
                self._run_basic_limiter_cpp(results)
            else:
                self._run_basic_limiter(results)
        
        if self.filter == 'minmax':
            self._run_minmax_filter(results)
        
        self.phi.vector().set_local(results, self.local_indices_dg1)
        self.phi.vector().apply('insert')
    
    def _run_basic_limiter_cpp(self, results):
        num_cells_all = self.mesh.num_cells()
        tdim = self.mesh.topology().dim()
        num_cells_owned = self.mesh.topology().ghost_offset(tdim)
        exceedances = self.excedance.vector().get_local()
        self.cpp_mod.slope_limiter_basic_cg1(self.num_neighbours,
                                             num_cells_all,
                                             num_cells_owned,
                                             self.max_neighbours,
                                             self.flat_neighbours,
                                             self.flat_cell_dofs,
                                             self.flat_cell_dofs_dg0,
                                             self.flat_vertices,
                                             exceedances,
                                             results)
        self.excedance.vector().set_local(exceedances)
        self.excedance.vector().apply('insert')
    
    def _run_basic_limiter(self, results):
        """
        Perform basic slope limiting
        """
        exceedances = self.excedance.vector().get_local()
        
        # Get number of non-ghosts
        assert self.mesh.topology().dim() ==  2
        ncells = self.mesh.topology().ghost_offset(2)
        
        # Cell averages
        averages = []
        for ic in range(self.mesh.num_cells()):
            dofs = self.cell_dofs_DG1[ic]
            vals = [results[dof] for dof in dofs]
            averages.append(sum(vals)/3)
        
        # Modify vertex values
        onetwothree = range(3)
        for ic in range(ncells):
            vertices = self.vertices[ic]
            dofs = self.cell_dofs_DG1[ic]
            vals = [results[dof] for dof in dofs]
            avg = sum(vals)/3
            
            excedance = 0
            for ivertex in onetwothree:
                vtx = vertices[ivertex]
                nbs = self.neighbours[vtx]
                nb_vals = [averages[self.cell_dofs_DG0[nb]] for nb in nbs]
                
                # Find highest and lowest value in the connected cells
                lo, hi = 1e100, -1e100
                for cell_avg in nb_vals:
                    lo = min(lo, cell_avg)
                    hi = max(hi, cell_avg)
                
                vtx_val = vals[ivertex] 
                if vtx_val < lo:
                    vals[ivertex] = lo
                    ex = vtx_val - lo
                    if abs(excedance) < abs(ex):
                        excedance = ex
                elif vtx_val > hi:
                    vals[ivertex] = hi
                    ex = vtx_val - hi
                    if abs(excedance) < abs(ex):
                        excedance = ex
            
            exceedances[self.cell_dofs_DG0[ic]] = excedance
            if excedance == 0:
                continue
            
            # Modify the results to limit the slope
            
            # Find the new average and which vertices can be adjusted to obtain the correct average
            new_avg = sum(vals)/3
            eps = 0
            moddable = [0, 0, 0]
            if abs(avg - new_avg) > 1e-15:
                if new_avg > avg:
                    for ivertex in onetwothree:
                        if vals[ivertex] > avg:
                            moddable[ivertex] = 1
                else:
                    for ivertex in onetwothree:
                        if vals[ivertex] < avg:
                            moddable[ivertex] = 1
                
                # Get number of vertex values that can be modified and catch
                # possible floating point problems with the above comparisons
                nmod = sum(moddable)
                
                if nmod == 0:
                    assert abs(excedance) < 1e-14, 'Nmod=0 with exceedance=%r' % excedance
                else:
                    eps = (avg - new_avg)*3/nmod
            
            # Modify the vertices to obtain the correct average
            for iv in onetwothree:
                results[dofs[iv]] = vals[iv] + eps*moddable[iv]
        
        self.excedance.vector().set_local(exceedances)
        self.excedance.vector().apply('insert')
    
    def _run_minmax_filter(self, results):
        """
        Make sure the DOF values are inside the min and max values from the
        first time step (enforce 
        """
        if self._filter_cache is None:
            mi = df.MPI.min(df.mpi_comm_world(), float(results.min()))
            ma = df.MPI.max(df.mpi_comm_world(), float(results.max()))
            self._filter_cache = mi, ma
            return
        minval, maxval = self._filter_cache
        numpy.clip(results, minval, maxval, results)


class DoNothingSlopeLimiter(object):
    def __init__(self, *argv, **kwargs):
        pass
    
    def run(self):
        pass


def SlopeLimiter(simulation, phi_name, phi):
    method = simulation.input.get_value('slope_limiter/%s/method', LIMITER, 'string')
    filter_method = simulation.input.get_value('slope_limiter/%s/filter', FILTER, 'string')
    use_cpp = simulation.input.get_value('slope_limiter/%s/use_cpp', USE_CPP, 'boolean')
    
    verify_key('slope limiter', method, _METHODS, '%s transport' % phi_name)
    
    simulation.log.info('    Using slope limiter %s with filter %s for %s' % (method, filter, phi_name))
    limiter = _METHODS[method](phi, filter_method, use_cpp)
    return limiter


_METHODS = {}
_METHODS['none'] = DoNothingSlopeLimiter
_METHODS['basic'] = SlopeLimiterBasic
