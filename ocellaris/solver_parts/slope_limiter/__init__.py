import numpy
import dolfin
from ocellaris.utils import ocellaris_error
from ocellaris.solver_parts import get_dof_region_marks


LIMITER = 'none'
FILTER = 'nofilter'
USE_CPP = True
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


@register_slope_limiter('None')
class DoNothingSlopeLimiter(SlopeLimiterBase):
    description = 'No slope limiter'
    
    def __init__(self, *argv, **kwargs):
        self.additional_plot_funcs = []
    
    def run(self):
        pass


def SlopeLimiter(simulation, phi_name, phi, default_limiter=LIMITER, default_filter=FILTER, default_use_cpp=USE_CPP):
    """
    Return a slope limiter based on the user provided input or the default
    values if no input is provided by the user
    """
    # Get user provided input (or default values)
    inp = simulation.input.get_value('slope_limiter/%s' % phi_name, {}, 'Input')
    method = inp.get_value('method', default_limiter, 'string')
    filter_method = inp.get_value('filter', default_filter, 'string')
    use_cpp = inp.get_value('use_cpp', default_use_cpp, 'bool')
    plot_exceedance = inp.get_value('plot', False, 'bool')  
    
    # Get the region markers
    V = phi.function_space()
    dof_region_marks = get_dof_region_marks(simulation, V)
    boundary_condition = numpy.zeros(V.dim(), numpy.intc)
    for dof in dof_region_marks:
        boundary_condition[dof] = 1
    
    # Construct the limiter
    simulation.log.info('    Using slope limiter %s with filter %s for %s' % (method, filter_method, phi_name))
    limiter_class = get_slope_limiter(method)
    limiter = limiter_class(phi_name, phi, boundary_condition, filter_method, use_cpp)
    
    if plot_exceedance:
        for func in limiter.additional_plot_funcs:
            simulation.io.add_extra_output_function(func)
    
    return limiter


from . import naive_nodal
from . import hierarchal_taylor
