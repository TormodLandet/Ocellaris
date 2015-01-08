import functools
import numpy
import dolfin
from dgvof.utils import report_error, GradientReconstructor

_CONVECTION_SCHEMES = {}

def add_convection_scheme(name, scheme_class):
    """
    Register a convection scheme
    """
    _CONVECTION_SCHEMES[name] = scheme_class

def register_convection_scheme(name):
    """
    A class decorator to register convection schemes
    """
    def register(scheme_class):
        add_convection_scheme(name, scheme_class)
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
    
    def __init__(self, alpha_function):
        """
        The given function space is for the function you
        will be convected. The convection scheme itself
        uses a Crouzeix-Raviart 1st order element to 
        represent the blending function 
        
        The blending function is a downstream blending
        factor (0=upstream, 1=downstream)
        
        The alpha function is the scalar function to be
        advected
        """ 
        self.alpha_function = alpha_function
        self.alpha_function_space = alpha_function.function_space()
        self.mesh = self.alpha_function_space.mesh()
        
        # Function space for the convection blending function
        self.function_space = dolfin.FunctionSpace(self.mesh, "CR", 1)
        self.blending_function = dolfin.Function(self.function_space)
        
        # Dofmaps
        self.alpha_dofmap = self.alpha_function_space.dofmap().dofs()
        self.dofmap = self.function_space.dofmap().dofs()
        
        # Connectivity from face to edge
        self.mesh.init(2, 1)
        self.con21 = self.mesh.topology()(2, 1)

        # Connectivity from edge to face
        self.mesh.init(1, 2)
        self.con12 = self.mesh.topology()(1, 2)
        
        # Topological dimension
        self.ndim = self.function_space.cell().topological_dimension()
        
        # Mesh size
        self.ncells = self.mesh.num_cells()
        self.nfacets = self.mesh.num_facets()
        
        # Cell centroids
        self.centroids = numpy.zeros((self.ncells, self.ndim), float)
        for cell in dolfin.cells(self.mesh):
            mp = cell.midpoint()
            if self.ndim == 2:
                self.centroids[cell.index()] = (mp.x(), mp.y())
            else:
                self.centroids[cell.index()] = (mp.x(), mp.y(), mp.z())

        self.facet_centroids = numpy.zeros((self.nfacets, self.ndim), float)
        for facet in dolfin.facets(self.mesh):
            mp = facet.midpoint()
            self.facet_centroids[facet.index()] = (mp.x(), mp.y())

        # To allow comparison
        self.force_upwind = False
        
        # For gradient reconstruction
        self.gradient_reconstructor = GradientReconstructor(self.alpha_function, self.alpha_dofmap, self.centroids)

    def update(self, t, dt, velocity):
        raise NotImplementedError()

from . import hric
