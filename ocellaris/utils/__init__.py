from .error_handling import OcellarisError, ocellaris_error, verify_key
from .interactive_console import interactive_console_hook, run_debug_console
from .timer import timeit, log_timings
from .code_runner import RunnablePythonString, CodedExpression
from .cpp_expression import OcellarisCppExpression, ocellaris_interpolate
from .gradient_reconstruction import GradientReconstructor
from .dofmap import facet_dofmap, get_dof_neighbours
from .linear_solvers import linear_solver_from_input, condition_number
from .trace_projection import convert_to_dgt
from .mpi import get_root_value, gather_lines_on_root
from .taylor_basis import lagrange_to_taylor, taylor_to_lagrange
from .small_helpers import create_vector_functions, shift_fields, velocity_change