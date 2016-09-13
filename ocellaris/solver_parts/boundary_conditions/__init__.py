import dolfin
from ocellaris.utils import ocellaris_error, RunnablePythonString


_BOUNDARY_CONDITIONS = {}


def add_boundary_condition(name, boundary_condition_class):
    """
    Register a boundary condition
    """
    _BOUNDARY_CONDITIONS[name] = boundary_condition_class


def register_boundary_condition(name):
    """
    A class decorator to register boundary conditions
    """
    def register(boundary_condition_class):
        add_boundary_condition(name, boundary_condition_class)
        return boundary_condition_class
    return register


def get_boundary_condition(name):
    """
    Return a boundary condition by name
    """
    try:
        return _BOUNDARY_CONDITIONS[name]
    except KeyError:
        ocellaris_error('Boundary condition "%s" not found' % name,
                        'Available boundary conditions:\n' +
                        '\n'.join('  %-20s - %s' % (n, s.description)
                                  for n, s in sorted(_BOUNDARY_CONDITIONS.items())))
        raise


class BoundaryCondition(object):
    description = 'No description available'
    
    def func(self):
        """
        Returns the value at the boundary for Dirichlet boundary conditions
        and the normal derivative at the boundaru for Neumann bcs.
        """
        raise NotImplementedError()
    
    def ds(self):
        """
        Returns the ds measure of the part of the boundary which this boundary
        condition applies to
        """
        raise NotImplementedError()


from .boundary_region import BoundaryRegion
from .dof_marker import get_dof_region_marks

from . import dirichlet
from . import neumann
from . import wall
from . import outlet
