import numpy
import dolfin
from ocellaris.utils import ocellaris_error
from ocellaris.solver_parts import get_dof_region_marks, mark_cell_layers
from .limiter_bcs import SlopeLimiterBoundaryConditions


DEFAULT_LIMITER = 'None'
DEFAULT_FILTER = 'nofilter'
DEFAULT_USE_CPP = True
_SLOPE_LIMITERS = {}


def add_slope_limiter(name, slope_limiter_class):
    """
    Register a slope limiter
    """
    _SLOPE_LIMITERS[name] = slope_limiter_class


def register_slope_limiter(name):
    """
    A class decorator to register slope limiters
    """
    def register(slope_limiter_class):
        add_slope_limiter(name, slope_limiter_class)
        slope_limiter_class.limiter_method = name 
        return slope_limiter_class
    return register


def get_slope_limiter(name):
    """
    Return a slope limiter by name
    """
    try:
        return _SLOPE_LIMITERS[name]
    except KeyError:
        ocellaris_error('Slope limiter "%s" not found' % name,
                        'Available slope limiters:\n' +
                        '\n'.join('  %-20s - %s' % (n, s.description)
                                  for n, s in sorted(_SLOPE_LIMITERS.items())))
        raise


class SlopeLimiterBase(object):
    description = 'No description available'
    
    active = True
    phi_old = None
    global_bounds = -1e100, 1e100
    has_global_bounds = False
    enforce_global_bounds = False
    
    def set_phi_old(self, phi_old):
        """
        Some limiters use an old (allready limited) field as an extra
        information about in the limiting process
        """
        self.phi_old = phi_old
    
    def set_global_bounds(self, phi=None, lo=None, hi=None):
        """
        Set the initial field which we use this to compute the global
        bounds which must be observed everywhere for all time 
        """
        if self.enforce_global_bounds:
            if phi is None:
                assert isinstance(lo, (int, float))
                assert isinstance(hi, (int, float))
            else:    
                lo = phi.vector().min()
                hi = phi.vector().max()
            self.global_bounds = (lo, hi)
            self.has_global_bounds = True
        else:
            self.global_bounds = SlopeLimiterBase.global_bounds
            self.has_global_bounds = False
        return self.global_bounds


@register_slope_limiter('None')
class DoNothingSlopeLimiter(SlopeLimiterBase):
    description = 'No slope limiter'
    active = False
    
    def __init__(self, *argv, **kwargs):
        self.additional_plot_funcs = []
    
    def run(self):
        pass


@register_slope_limiter('OnlyBound')
class OnlyBoundSlopeLimiter(SlopeLimiterBase):
    description = 'Bounding limiter'
    active = False
    
    def __init__(self, phi_name, phi, enforce_bounds=False, **kwargs):
        self.phi_name = phi_name
        self.phi = phi
        self.enforce_global_bounds = enforce_bounds
        self.additional_plot_funcs = []
    
    def run(self):
        if self.has_global_bounds:
            lo, hi = self.global_bounds
            vals = self.phi.vector().get_local()
            numpy.clip(vals, lo, hi, vals)
            self.phi.vector().set_local(vals)
            self.phi.vector().apply('insert')


def SlopeLimiter(simulation, phi_name, phi, output_name=None, method=None):
    """
    Return a slope limiter based on the user provided input or the default
    values if no input is provided by the user
    """
    # Get user provided input (or default values)
    inp = simulation.input.get_value('slope_limiter/%s' % phi_name, {}, 'Input')
    if method is None:
        method = inp.get_value('method', DEFAULT_LIMITER, 'string')
    use_cpp = inp.get_value('use_cpp', DEFAULT_USE_CPP, 'bool')
    plot_exceedance = inp.get_value('plot', False, 'bool')
    skip_boundary = inp.get_value('skip_boundary', False, 'bool')
    enforce_bounds = inp.get_value('enforce_bounds', False, 'bool')
    enforce_bcs = inp.get_value('enforce_bcs', True, 'bool')
    name = phi_name if output_name is None else output_name
    
    # Find degree
    V = phi.function_space()
    degree = V.ufl_element().degree()
    
    # Skip limiting if degree is 0
    if degree == 0:
        if 'method' in inp and method != 'OnlyBound':
            simulation.log.info('    Switching to slope limiter OnlyBound for field '
                                '%s (due to degree == 0)' % name)
        return OnlyBoundSlopeLimiter(phi_name=phi_name, phi=phi, enforce_bounds=enforce_bounds)
    
    # Get boundary region marks and get the helper class used to limit along the boundaries
    drm = get_dof_region_marks(simulation, V)
    bcs = SlopeLimiterBoundaryConditions(simulation, phi_name, drm, V)
    
    if skip_boundary:
        # Mark boundary cells for skipping (cells containing dofs with region marks)
        skip_cells = mark_cell_layers(simulation, V, layers=0, dof_region_marks=drm)
    else:
        skip_cells = ()
        bcs.activate()
    
    # Construct the limiter
    simulation.log.info('    Using slope limiter %s for field %s' % (method, name))
    simulation.log.info('        Enforcing BCs: %r' % enforce_bcs)
    simulation.log.info('        Skip boundary: %r' % skip_boundary)
    simulation.log.info('        Enforce global bounds: %r' % enforce_bounds)
    limiter_class = get_slope_limiter(method)
    limiter = limiter_class(phi_name=phi_name, phi=phi, skip_cells=skip_cells,
                            boundary_conditions=bcs, output_name=output_name,
                            use_cpp=use_cpp, enforce_bounds=enforce_bounds,
                            enforce_bcs=enforce_bcs)
    
    if plot_exceedance:
        for func in limiter.additional_plot_funcs:
            simulation.io.add_extra_output_function(func)
    
    return limiter


from ocellaris.cpp import load_module
LocalMaximaMeasurer = load_module('measure_local_maxima').LocalMaximaMeasurer

from . import naive_nodal
from . import hierarchical_taylor
