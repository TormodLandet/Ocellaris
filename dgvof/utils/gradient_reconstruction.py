import dolfin
import numpy

class GradientReconstructor(object):
    def __init__(self, simulation, alpha_func, alpha_dofmap,
                 use_vertex_neighbours=True, compute_facet_gradient=False):
        """
        Reconstructor for the gradient in each cell. Assumes DG0
        
        See for example "An Introduction to Computational Fluid Dynamics -
        The Finite Volume Method" by Versteeg & Malalasekera (2007), 
        specifically equation 11.36 on page 322 for details on the method
        """
        self.simulation = simulation
        self.alpha_function = alpha_func
        self.alpha_dofmap = alpha_dofmap
        self.mesh = alpha_func.function_space().mesh()
        self.use_vertex_neighbours = use_vertex_neighbours
        self.reconstruction_initialized = False
        
        self.compute_facet_gradient = compute_facet_gradient
        if self.compute_facet_gradient:
            # The gradient of the alpha function at the faces
            V = dolfin.VectorFunctionSpace(self.mesh, "CR", 1)
            self.facet_gradient = dolfin.Function(V)
            
            # The dofmap of this vector function
            self.facet_gradient_dofmap = []
            ndim  = V.cell().topological_dimension()
            for i in xrange(ndim):
                self.facet_gradient_dofmap.append(V.sub(i).dofmap().dofs())
    
    def initialize(self):
        """
        Precompute least squares matrices
        """
        V = self.alpha_function.function_space()
        Vvec = dolfin.VectorFunctionSpace(self.mesh, 'DG', 0)
        ndim = V.cell().topological_dimension()
        ncells = self.mesh.num_cells()
        
        self.gradient = dolfin.Function(Vvec)
        self.gradient_dofmap0 = Vvec.sub(0).dofmap().dofs()
        self.gradient_dofmap1 = Vvec.sub(1).dofmap().dofs()
        self.neighbour_minval = dolfin.Function(V)
        self.neighbour_maxval = dolfin.Function(V)
        self.neighbours = [None]*ncells
        self.lstsq_matrices = [None]*ncells
        self.lstsq_inv_matrices = numpy.zeros((ncells, ndim, ndim), float)
        
        cell_info = self.simulation.data['cell_info']
        conFC = self.simulation.data['connectivity_FC']
        conCF = self.simulation.data['connectivity_CF']
        conCC = self.simulation.data['connectivity_CC']
        
        for cell in dolfin.cells(self.mesh):
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
            dix = self.alpha_dofmap[idx]
            neighbours = self.neighbours[idx]
            
            # Get the matrices
            AT = self.lstsq_matrices[idx]
            ATAI = self.lstsq_inv_matrices[idx]
            a0  = a_cell_vec[self.alpha_dofmap[idx]]
            b = [(a_cell_vec[self.alpha_dofmap[ni]] - a0) for ni in neighbours]
            b = numpy.array(b, float)
            
            # Store min and max values which can be used to enforce convective boundedness
            self.neighbour_minval.vector()[dix] = b.min()
            self.neighbour_maxval.vector()[dix] = b.max()
            
            # Calculate the and store the gradient
            g = numpy.dot(ATAI, numpy.dot(AT, b))
            self.gradient.vector()[self.gradient_dofmap0[idx]] = g[0]
            self.gradient.vector()[self.gradient_dofmap1[idx]] = g[1] 
        
        # Calculate the gradient of the alpha field on the facets
#        if self.compute_facet_gradient:
#            V = self.facet_gradient.function_space()
#            self.facet_gradient.assign(dolfin.project(self.gradient, V=V))
#            
#             facet_gradient_vec = self.facet_gradient.vector()
#             fdofs = self.facet_gradient_dofmap
#             for facet in dolfin.facets(self.mesh):
#                 fidx = facet.index()
#                 
#                 # Find the local cells (the two cells sharing this face)
#                 connected_cells = self.con12(fidx)
#     
#                 if len(connected_cells) == 1:
#                     # Indices of the connected cell
#                     idx0, = connected_cells
#                     
#                     # Aproximate gradient at the interface
#                     g = self.gradient[idx0]
#                 
#                 else:
#                     # Indices of the two local cells
#                     idx0, idx1 = connected_cells
#                     
#                     # Aproximate gradient at the interface
#                     grad = 0.5*(self.gradient[idx0] + self.gradient[idx1])
#                 
#                 for i, g in enumerate(grad):
#                     facet_gradient_vec[fdofs[i][fidx]] = g
        