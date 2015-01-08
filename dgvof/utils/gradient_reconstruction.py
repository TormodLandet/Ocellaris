import dolfin
import numpy

class GradientReconstructor(object):
    def __init__(self, alpha_func, alpha_dofmap, centroids, use_vertex_neighbours=True):
        """
        Reconstructor for the gradient in each cell. Assumes DG0
        
        See for example "An Introduction to Computational Fluid Dynamics -
        The Finite Volume Method" by Versteeg & Malalasekera (2007), 
        specifically equation 11.36 on page 322 for details on the method
        """
        self.alpha_function = alpha_func
        self.alpha_dofmap = alpha_dofmap
        self.centroids = centroids
        self.use_vertex_neighbours = use_vertex_neighbours
        self.reconstruction_initialized = False
        
    def initialize(self):
        """
        Precompute least squares matrices
        """
        v = self.alpha_function.function_space()
        self.mesh = v.mesh()
        ndim = v.cell().topological_dimension()
        ncells = self.mesh.num_cells()
        
        self.gradient = numpy.zeros((ncells, ndim), float)
        self.neighbour_minval = numpy.zeros(ncells, float)
        self.neighbour_maxval = numpy.zeros(ncells, float)
        self.neighbours = [None]*ncells
        self.lstsq_matrices = [None]*ncells
        self.lstsq_inv_matrices = numpy.zeros((ncells, ndim, ndim), float)
        
        if self.use_vertex_neighbours:
            # Connectivity from face to face
            self.mesh.init(2, 2)
            con22 = self.mesh.topology()(2, 2)
        else:
            # Connectivity from face to edge
            self.mesh.init(2, 1)
            con21 = self.mesh.topology()(2, 1)
            
            # Connectivity from edge to face
            self.mesh.init(1, 2)
            con12 = self.mesh.topology()(1, 2)
        
        for cell in dolfin.cells(self.mesh):
            idx = cell.index()
            
            # Find neighbours
            if self.use_vertex_neighbours:
                # Find cells sharing one or more vertices
                neighbours = con22(idx)
            else:
                # Find cells sharing one or more facets
                neighbours = []
                facets = con21(idx)
                for fi in facets:
                    cells = con12(fi)
                    for cell in cells:
                        neighbours.extend([ci for ci in cells if ci != idx and ci not in neighbours])
            
            # Get the centroid of the cell neighbours
            nneigh = len(neighbours)
            A = numpy.zeros((nneigh, ndim), float)
            for j, ni in enumerate(neighbours):
                A[j] = self.centroids[ni] - self.centroids[idx]
            
            # Calculate the matrices needed for least squares gradient reconstruction
            AT = A.T
            ATA = numpy.dot(AT, A)
            self.neighbours[idx] = neighbours
            self.lstsq_matrices[idx] = AT
            self.lstsq_inv_matrices[idx] = numpy.linalg.inv(ATA)
            
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
        
        # Reconstruct the gradient for 
        a_cell_vec = self.alpha_function.vector()
        for cell in dolfin.cells(self.mesh):
            idx = cell.index()
            neighbours = self.neighbours[idx]
            
            # Get the matrices
            AT = self.lstsq_matrices[idx]
            ATAI = self.lstsq_inv_matrices[idx]
            a0  = a_cell_vec[self.alpha_dofmap[idx]]
            b = [(a_cell_vec[self.alpha_dofmap[ni]] - a0) for ni in neighbours]
            b = numpy.array(b, float)
            
            # Store min and max values which can be used to enforce convective boundedness
            self.neighbour_minval[idx] = b.min()
            self.neighbour_maxval[idx] = b.max()
            
            # Calculate the and store the gradient
            self.gradient[idx] = numpy.dot(ATAI, numpy.dot(AT, b))
