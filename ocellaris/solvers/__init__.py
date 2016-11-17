import dolfin
from ocellaris.utils import ocellaris_error


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
        ocellaris_error('Navier-Stokes solver "%s" not found' % name,
                        'Available solvers:\n' +
                        '\n'.join('  %-20s - %s' % (n, s.description) 
                                  for n, s in sorted(_SOLVERS.items())))
        raise


class Solver(object):
    description = 'No description available'
    

class BaseEquation(object):
    # Will be shadowed by object properties after first assemble
    tensor_lhs = None
    tensor_rhs = None
    
    def assemble_lhs(self):
        if self.tensor_lhs is None:
            self.tensor_lhs = dolfin.assemble(self.form_lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            self.tensor_rhs = dolfin.assemble(self.form_rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs


# Timestepping methods
BDF = 'BDF'
CRANK_NICOLSON = 'CN'


# Flux types
BLENDED = 'Blended'
UPWIND = 'Upwind'
LOCAL_LAX_FRIEDRICH = 'Local Lax-Friedrich'


# Velocity post-processing
BDM = 'BDM'


from . import analytical_solution
from . import coupled
from . import coupled_ldg
from . import fsvd
from . import ipcs
from . import simple
