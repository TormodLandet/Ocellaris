import functools
import numpy
import dolfin
from dgvof.utils import report_error

_CONVECTION_SCHEMES = {}

def register_convection_scheme(name, scheme_class):
    """
    Register a convection scheme
    """
    _CONVECTION_SCHEMES[name] = scheme_class

def convection_scheme(name):
    """
    A class decorator to register convection schemes
    """
    def register(scheme_class):
        register_convection_scheme(name, scheme_class)
        return scheme_class
    return register

def get_convection_scheme(name):
    """
    Return a convection scheme by name
    """
    try:
        return _CONVECTION_SCHEMES[name]
    except KeyError:
        report_error('Convection scheme "%s" not found' % name,
                     'Available convection schemes:\n' +
                     '\n'.join('  %-20s - %s' % (n, s.description) 
                               for n, s in sorted(_CONVECTION_SCHEMES.items())))
        raise

class ConvectionScheme(object):
    """
    A generic convection scheme, to be subclassed
    for actual usage
    """
    description = 'No description available'
    
    def __init__(self, mesh, alpha_function_space):
        """
        The given function space is for the function you
        will be convected. The convection scheme itself
        uses a Crouzeix-Raviart 1st order element to 
        represent the blending function 
        
        The blending function is a downstream blending
        factor (0=upstream, 1=downstream)
        """
        self.mesh = mesh
        
        # Function space for the function to be advected
        self.alpha_function_space = alpha_function_space
        
        # Function space for the convection blending function
        self.function_space = dolfin.FunctionSpace(mesh, "CR", 1)
        self.blending_function = dolfin.Function(self.function_space)
        
        # Dofmaps
        self.alpha_dofmap = self.alpha_function_space.dofmap().dofs()
        self.dofmap = self.function_space.dofmap().dofs()
        
        # Connectivity from face to edge
        mesh.init(2, 1)
        self.con21 = mesh.topology()(2, 1)

        # Connectivity from edge to face
        mesh.init(1, 2)
        self.con12 = mesh.topology()(1, 2)
        
        # Topological dimension
        self.ndim = self.function_space.cell().topological_dimension()
        
        # Mesh size
        self.ncells = mesh.num_cells()
        self.nfacets = mesh.num_facets()
        
        # Cell centroids
        self.centroids = numpy.zeros((self.ncells, self.ndim), float)
        for cell in dolfin.cells(mesh):
            mp = cell.midpoint()
            if self.ndim == 2:
                self.centroids[cell.index()] = (mp.x(), mp.y())
            else:
                self.centroids[cell.index()] = (mp.x(), mp.y(), mp.z())

        self.facet_centroids = numpy.zeros((self.nfacets, self.ndim), float)
        for facet in dolfin.facets(mesh):
            mp = facet.midpoint()
            self.facet_centroids[facet.index()] = (mp.x(), mp.y())

        # To allow comparison
        self.force_upwind = False
        
        # For gradient reconstruction
        self.gradient_reconstruction_initialized = False
        self.gradient = numpy.zeros((self.ncells, self.ndim), float)
        self.neighbour_maxval = numpy.zeros(self.ncells, float)
        self.neighbour_minval = numpy.zeros(self.ncells, float)

    def update(self, t, dt, alpha_function, velocity):
        raise NotImplementedError()
    
    def reconstruct_gradient(self, alpha_function):
        """
        Reconstruct the gradient in each cell. Assumes DG0
        
        TODO: handle boundary conditions for boundary cells,
              right now the boundary cell gradients are only
              influenced by the cell neighbours
        """
        # Initialize the least squares gradient reconstruction matrices
        # See Verstee & Malalasekera (2007) equation 11.36, p 322
        if not self.gradient_reconstruction_initialized:
            self.gradient_neighbours = [None]*self.ncells
            self.gradient_lstsq_matrices = [None]*self.ncells
            self.gradient_lstsq_inv_matrices = numpy.zeros((self.ncells, self.ndim, self.ndim), float) 
            for cell in dolfin.cells(self.mesh):
                idx = cell.index()
                
                # Find cell neighbours
                neighbours = []
                facets = self.con21(idx)
                for fi in facets:
                    cells = self.con12(fi)
                    for cell in cells:
                        neighbours.extend([ci for ci in cells if ci != idx and ci not in neighbours])
                
                # Get the centroid of the cell neighbours
                nneigh = len(neighbours)
                A = numpy.zeros((nneigh, self.ndim), float)
                for j, ni in enumerate(neighbours):
                    A[j] = self.centroids[ni] - self.centroids[idx]
                
                # Calculate the matrices needed for least squares gradient reconstruction
                AT = A.T
                ATA = numpy.dot(AT, A)
                self.gradient_neighbours[idx] = neighbours
                self.gradient_lstsq_matrices[idx] = AT
                self.gradient_lstsq_inv_matrices[idx] = numpy.linalg.inv(ATA)
                
            self.gradient_reconstruction_initialized = True
        
        # Reconstruct the gradient for 
        a_cell_vec = alpha_function.vector()
        for cell in dolfin.cells(self.mesh):
            idx = cell.index()
            neighbours = self.gradient_neighbours[idx]
            
            # Get the matrices
            AT = self.gradient_lstsq_matrices[idx]
            ATAI = self.gradient_lstsq_inv_matrices[idx]
            b = [(a_cell_vec[ni] - a_cell_vec[idx]) for ni in neighbours]
            b = numpy.array(b, float)
            
            # Store min and max values which can be used to enforce convective boundedness
            self.neighbour_maxval[idx] = b.max()
            self.neighbour_maxval[idx] = b.min()
            
            # Calculate the gradient
            self.gradient[idx] = numpy.dot(ATAI, numpy.dot(AT, b))

from . import hric
