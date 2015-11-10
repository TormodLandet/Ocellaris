import dolfin
import numpy
from ocellaris.cpp import load_module


class GradientReconstructor(object):
    def __init__(self, simulation, alpha_func, use_vertex_neighbours=True):
        """
        Reconstructor for the gradient in each cell.
        
        See for example "An Introduction to Computational Fluid Dynamics -
        The Finite Volume Method" by Versteeg & Malalasekera (2007), 
        specifically equation 11.36 on page 322 for details on the method
        """
        self.simulation = simulation
        self.alpha_function = alpha_func
        self.mesh = alpha_func.function_space().mesh()
        self.use_vertex_neighbours = use_vertex_neighbours
        self.reconstruction_initialized = False
        self.use_cpp = simulation.input.get_value('use_cpp_extensions', True, 'bool')
    
    def initialize(self):
        """
        Initialize the gradient function and dofmap
        
        For DG0 fields we must also precompute the least squares matrices
        needed to reconstruct gradients since the gradient of a constant
        is zero in each cell so a FEM projection will not work
        """
        V = self.alpha_function.function_space()
        cd = self.simulation.data['constrained_domain']
        Vvec = dolfin.VectorFunctionSpace(self.mesh, 'DG', 0, constrained_domain=cd)
        ndim = V.ufl_cell().topological_dimension()
        ncells = self.mesh.num_cells()
        
        # To be used by others accessing this class
        self.gradient = dolfin.Function(Vvec)
        self.gradient_dofmaps = [Vvec.sub(d).dofmap().dofs() for d in range(ndim)]
        
        if self.alpha_function.ufl_element().degree() > 0:
            # We do not need the rest of the precomputed data for
            # higher order functions
            return
        
        # Connectivity info needed in calculations
        cell_info = self.simulation.data['cell_info']
        conFC = self.simulation.data['connectivity_FC']
        conCF = self.simulation.data['connectivity_CF']
        conVC = self.simulation.data['connectivity_VC']
        conCV = self.simulation.data['connectivity_CV']
        
        if self.use_vertex_neighbours:
            # Find cells sharing one or more vertices
            con1 = conCV
            con2 = conVC
        else:
            # Find cells sharing one or more facets
            con1 = conCF
            con2 = conFC
        
        # Precompute connectivity and geometry matrices 
        everyones_neighbours = [None]*ncells
        lstsq_matrices = [None]*ncells
        self.lstsq_inv_matrices = numpy.zeros((ncells, ndim, ndim), float, order='C')
        
        for cell in dolfin.cells(self.mesh):
            idx = cell.index()
            
            # Find neighbours
            neighbours = []
            facets_or_vertices = con1(idx)
            for ifv in facets_or_vertices:
                cell_neighbours = con2(ifv)
                neighbours.extend([ci for ci in cell_neighbours if ci != idx and ci not in neighbours])
            
            # Get the centroid of the cell neighbours
            nneigh = len(neighbours)
            A = numpy.zeros((nneigh, ndim), float)
            mp0 = cell_info[idx].midpoint
            for j, ni in enumerate(neighbours):
                mpJ = cell_info[ni].midpoint
                A[j] = mpJ - mp0 
            
            # Calculate the matrices needed for least squares gradient reconstruction
            AT = A.T
            ATA = numpy.dot(AT, A)
            everyones_neighbours[idx] = neighbours
            lstsq_matrices[idx] = AT
            self.lstsq_inv_matrices[idx] = numpy.linalg.inv(ATA)
        
        # Turn the lists into numpy arrays for ease of communication with C++
        N = len(everyones_neighbours)
        self.num_neighbours = numpy.array([len(nbs) for nbs in everyones_neighbours], dtype='i', order='C')
        NBmax = self.num_neighbours.max()
        self.neighbours = numpy.zeros((N, NBmax), dtype='i', order='C')
        self.lstsq_matrices = numpy.zeros((N, ndim, NBmax), float, order='C')
        for i in xrange(N):
            Nnb = self.num_neighbours[i]
            self.neighbours[i,:Nnb] = everyones_neighbours[i]
            self.lstsq_matrices[i,:,:Nnb] = lstsq_matrices[i] 
        
        # Instant only allows one dimensional arrays
        self.lstsq_matrices = self.lstsq_matrices.reshape(-1, order='C')
        self.lstsq_inv_matrices = self.lstsq_inv_matrices.reshape(-1, order='C')
        self.neighbours = self.neighbours.reshape(-1, order='C')
        self.max_neighbours = NBmax
        
        self.reconstruction_initialized = True
    
    def reconstruct(self):
        """
        Reconstruct the gradient in each cell center
        
        TODO: handle boundary conditions for boundary cells,
              right now the boundary cell gradients are only
              influenced by the cell neighbours
        """
        # Initialize the least squares gradient reconstruction matrices
        # needed to calculate the gradient of a DG0 field
        if not self.reconstruction_initialized:
            self.initialize()
        
        if self.alpha_function.ufl_element().degree() > 0:
            # We use projection for higher degrees
            V = self.gradient.function_space()
            dolfin.project(dolfin.nabla_grad(self.alpha_function), V, function=self.gradient)
        else:
            if not self.use_cpp:
                # Pure Python version
                reconstructor = _reconstruct_gradient 
            else:
                # Faster C++ version
                cpp_gradient_reconstruction = load_module('gradient_reconstruction')
                reconstructor = cpp_gradient_reconstruction.reconstruct_gradient
            
            # Run the gradient reconstruction
            reconstructor(self.alpha_function,
                          self.num_neighbours,
                          self.max_neighbours,
                          self.neighbours, 
                          self.lstsq_matrices,
                          self.lstsq_inv_matrices,
                          self.gradient)


def _reconstruct_gradient(alpha_function, num_neighbours, max_neighbours, neighbours, lstsq_matrices, lstsq_inv_matrices, gradient):
    """
    Reconstruct the gradient, Python version of the code
    
    This function used to have a more Pythonyc implementation
    that was most likely also faster. See old commits for that
    code. This code is here to verify the C++ version that is
    much faster than this (and the old Pythonic version)
    """
    a_cell_vec = alpha_function.vector()
    mesh = alpha_function.function_space().mesh()
    
    V = alpha_function.function_space()
    alpha_dofmap = V.dofmap().dofs()
    Vvec = gradient.function_space()
    gradient_dofmap0 = Vvec.sub(0).dofmap().dofs()
    gradient_dofmap1 = Vvec.sub(1).dofmap().dofs()
    
    np_gradient = gradient.vector().get_local()

    # Reshape arrays. The C++ version needs flatt arrays
    # (limitation in Instant/Dolfin) and we have the same
    # interface for both versions of the code
    ncells = len(num_neighbours)
    ndim = mesh.topology().dim()
    neighbours = neighbours.reshape((ncells, max_neighbours))
    lstsq_matrices = lstsq_matrices.reshape((ncells, ndim, max_neighbours))
    lstsq_inv_matrices = lstsq_inv_matrices.reshape((ncells, ndim, ndim))
    
    for i, cell in enumerate(dolfin.cells(mesh)):
        idx = cell.index()
        dix = alpha_dofmap[idx]
        Nnbs = num_neighbours[i]
        nbs = neighbours[i,:Nnbs]
        
        # Get the matrices
        AT = lstsq_matrices[i,:,:Nnbs]
        ATAI = lstsq_inv_matrices[i]
        a0  = a_cell_vec[dix]
        b = [(a_cell_vec[alpha_dofmap[ni]] - a0) for ni in nbs]
        b = numpy.array(b, float)
        
        # Calculate the and store the gradient
        g = numpy.dot(ATAI, numpy.dot(AT, b))
        np_gradient[gradient_dofmap0[idx]] = g[0]
        np_gradient[gradient_dofmap1[idx]] = g[1]
    
    gradient.vector().set_local(np_gradient)
    gradient.vector().apply('insert')
