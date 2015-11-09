from ocellaris.utils import report_error

_SOLVERS = {}

def add_solver(name, solver_class):
    """
    Register a Navier-Stokes solver
    """
    _SOLVERS[name] = solver_class

def register_solver(name):
    """
    A class decorator to register Navier-Stokes solvers
    """
    def register(solver_class):
        add_solver(name, solver_class)
        return solver_class
    return register

def get_solver(name):
    """
    Return a Navier-Stokes solver by name
    """
    try:
        return _SOLVERS[name]
    except KeyError:
        report_error('Navier-Stokes solver "%s" not found' % name,
                     'Available solvers:\n' +
                     '\n'.join('  %-20s - %s' % (n, s.description) 
                               for n, s in sorted(_SOLVERS.items())),
                     stop=True)
        raise

class Solver(object):
    description = 'No description available'

# Timestepping methods
BDF = 'BDF'
CRANK_NICOLSON = 'CN'

# Flux types
BLENDED = 'Blended'
UPWIND = 'Upwind'
LOCAL_LAX_FRIEDRICH = 'Local Lax-Friedrich'

# Velocity post-processing
BDM = 'BDM'

from . import ipcs
from . import coupled
