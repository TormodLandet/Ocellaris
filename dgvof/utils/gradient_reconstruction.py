import dolfin
import numpy
from dgvof.cpp import load_module

class GradientReconstructor(object):
    def __init__(self, simulation, alpha_func, use_vertex_neighbours=True):
        """
        Reconstructor for the gradient in each cell. Assumes DG0
        
        See for example "An Introduction to Computational Fluid Dynamics -
        The Finite Volume Method" by Versteeg & Malalasekera (2007), 
        specifically equation 11.36 on page 322 for details on the method
        """
        self.simulation = simulation
        self.alpha_function = alpha_func
        self.mesh = alpha_func.function_space().mesh()
        self.use_vertex_neighbours = use_vertex_neighbours
        self.reconstruction_initialized = False
    
    def initialize(self):
        """
        Precompute least squares matrices
        """
        V = self.alpha_function.function_space()
        Vvec = dolfin.VectorFunctionSpace(self.mesh, 'DG', 0)
        ndim = V.cell().topological_dimension()
        ncells = self.mesh.num_cells()
        
        # To be used by others accessing this class
        self.gradient = dolfin.Function(Vvec)
        self.gradient_dofmap0 = Vvec.sub(0).dofmap().dofs()
        self.gradient_dofmap1 = Vvec.sub(1).dofmap().dofs()
        self.neighbour_minval = dolfin.Function(V)
        self.neighbour_maxval = dolfin.Function(V)
        
        # Connectivity info needed in calculations
        cell_info = self.simulation.data['cell_info']
        conFC = self.simulation.data['connectivity_FC']
        conCF = self.simulation.data['connectivity_CF']
        conCC = self.simulation.data['connectivity_CC']
         
        # Precompute connectivity and geometry matrices 
        everyones_neighbours = [None]*ncells
        lstsq_matrices = [None]*ncells
        self.lstsq_inv_matrices = numpy.zeros((ncells, ndim, ndim), float, order='C')
        
        for i, cell in enumerate(dolfin.cells(self.mesh)):
            idx = cell.index()
            
            # Find neighbours
            if self.use_vertex_neighbours:
                # Find cells sharing one or more vertices
                neighbours = conCC(idx)
            else:
                # Find cells sharing one or more facets
                neighbours = []
                facets = conCF(idx)
                for fi in facets:
                    cells = conFC(fi)
                    for cell in cells:
                        neighbours.extend([ci for ci in cells if ci != idx and ci not in neighbours])
            
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
            everyones_neighbours[i] = neighbours
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
              
        The following properties are accessible after reconstruction:
         - gradient - array of (ncells, 2) values
         - neighbour_minval - array of (ncells,) values
         - neighbour_maxval - array of (ncells,) values
        """
        # Initialize the least squares gradient reconstruction matrices
        if not self.reconstruction_initialized:
            self.initialize()
        
        if True:
            reconstructor = _reconstruct_gradient 
        else:
            cpp_gradient_reconstruction = load_module('gradient_reconstruction')
            reconstructor = cpp_gradient_reconstruction.reconstruct_gradient
        
        # Run the gradient reconstruction
        reconstructor(self.alpha_function,
                      self.num_neighbours,
                      self.max_neighbours,
                      self.neighbours, 
                      self.lstsq_matrices,
                      self.lstsq_inv_matrices,
                      self.neighbour_minval,
                      self.neighbour_maxval,
                      self.gradient)

def _reconstruct_gradient(alpha_function, num_neighbours, max_neighbours, neighbours, lstsq_matrices, lstsq_inv_matrices, neighbour_minval, neighbour_maxval, gradient):
    """
    Reconstruct the gradient, Python version of the code
    """
    a_cell_vec = alpha_function.vector()
    mesh = alpha_function.function_space().mesh()
    
    V = alpha_function.function_space()
    alpha_dofmap = V.dofmap().dofs()
    Vvec = gradient.function_space()
    gradient_dofmap0 = Vvec.sub(0).dofmap().dofs()
    gradient_dofmap1 = Vvec.sub(1).dofmap().dofs()
    
    np_gradient = gradient.vector().array()
    #np_minvals = neighbour_minval.vector().array()
    #np_maxvals = neighbour_maxval.vector().array()
    
    # Reshape arrays
    ncells = len(num_neighbours)
    ndim = mesh.topology().dim()
    neighbours = neighbours.reshape((ncells, max_neighbours))
    lstsq_matrices = lstsq_matrices.reshape((ncells, ndim, max_neighbours))
    lstsq_inv_matrices = lstsq_inv_matrices.reshape((ncells, ndim, ndim))
    
    for i, cell in enumerate(dolfin.cells(mesh)):
        idx = cell.index()
        #dix = alpha_dofmap[idx]
        Nnbs = num_neighbours[i]
        nbs = neighbours[i,:Nnbs]
        
        # Get the matrices
        AT = lstsq_matrices[i,:,:Nnbs]
        ATAI = lstsq_inv_matrices[i]
        a0  = a_cell_vec[alpha_dofmap[idx]]
        b = [(a_cell_vec[alpha_dofmap[ni]] - a0) for ni in nbs]
        b = numpy.array(b, float)
        
        # Store min and max values which can be used to enforce convective boundedness
        #np_minvals[dix] = b.min()
        #np_maxvals[dix] = b.max()
        
        # Calculate the and store the gradient
        g = numpy.dot(ATAI, numpy.dot(AT, b))
        np_gradient[gradient_dofmap0[idx]] = g[0]
        np_gradient[gradient_dofmap1[idx]] = g[1]
    
    gradient.vector()[:] = np_gradient
    #neighbour_minval.vector()[:] = np_minvals
    #neighbour_maxval.vector()[:] = np_maxvals 
        