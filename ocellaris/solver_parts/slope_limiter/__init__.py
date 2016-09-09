import dolfin
from ocellaris.utils import ocellaris_error


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


@register_slope_limiter('none')
class DoNothingSlopeLimiter(SlopeLimiterBase):
    description = 'No slope limiter'
    
    def __init__(self, *argv, **kwargs):
        pass
    
    def run(self):
        pass


def SlopeLimiter(simulation, phi_name, phi, default_limiter=LIMITER, default_filter=FILTER, default_use_cpp=USE_CPP):
    inp = simulation.input.get_value('slope_limiter/%s' % phi_name, {}, 'Input')
    method = inp.get_value('method', default_limiter, 'string')
    filter_method = inp.get_value('filter', default_filter, 'string')
    use_cpp = inp.get_value('use_cpp', default_use_cpp, 'boolean')
    
    simulation.log.info('    Using slope limiter %s with filter %s for %s' % (method, filter_method, phi_name))
    limiter_class = get_slope_limiter(method)
    limiter = limiter_class(phi, filter_method, use_cpp)
    return limiter


from . import basic_slope_limiter
