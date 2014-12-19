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
        
        # Cell centroids
        self.centroids = numpy.zeros((mesh.size(2), 2), float)
        for cell in dolfin.cells(mesh):
            mp = cell.midpoint()
            self.centroids[cell.index()] = (mp.x(), mp.y())

        self.facet_centroids = numpy.zeros((mesh.size(1), 2), float)
        for facet in dolfin.facets(mesh):
            mp = facet.midpoint()
            self.facet_centroids[facet.index()] = (mp.x(), mp.y())

        # To allow comparison
        self.force_upwind = False

    def update(self, t, dt, alpha_function, velocity):
        raise NotImplementedError()

from . import hric
