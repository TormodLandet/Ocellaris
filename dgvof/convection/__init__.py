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
    
    def __init__(self, simulation, func_name):
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
        self.simulation = simulation
        self.alpha_function = simulation.data[func_name]
        self.alpha_function_space = self.alpha_function.function_space()
        self.mesh = self.alpha_function_space.mesh()
        
        # Function space for the convection blending function
        self.function_space = dolfin.FunctionSpace(self.mesh, "CR", 1)
        self.blending_function = dolfin.Function(self.function_space)
        
        # Dofmaps
        self.alpha_dofmap = self.alpha_function_space.dofmap().dofs()
        self.dofmap = self.function_space.dofmap().dofs()
        
        # Mesh size
        self.ncells = self.mesh.num_cells()
        self.nfacets = self.mesh.num_facets()
        
        # For gradient reconstruction
        compute_facet_gradient = simulation.input.get('convection', {}).get(func_name, {}).get('compute_facet_gradient', False)
        self.gradient_reconstructor = GradientReconstructor(simulation, self.alpha_function, self.alpha_dofmap,
                                                            compute_facet_gradient=compute_facet_gradient)

    def update(self, t, dt, velocity):
        raise NotImplementedError()

from . import upwind
from . import cicsam
from . import hric
