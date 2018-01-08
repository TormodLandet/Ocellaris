from .error_handling import OcellarisError, ocellaris_error, verify_key
from .alarm import AlarmTimeout
from .interactive_console import interactive_console_hook, run_debug_console
from .timer import timeit, log_timings
from .code_runner import RunnablePythonString, CodedExpression
from .cpp_expression import OcellarisCppExpression, ocellaris_interpolate
from .dofmap import cell_dofmap, facet_dofmap, get_dof_neighbours
from .linear_solvers import (linear_solver_from_input, condition_number,
                             create_block_matrix, matmul)
from .mpi import get_root_value, sync_arrays, gather_lines_on_root
from .taylor_basis import lagrange_to_taylor, taylor_to_lagrange
from .small_helpers import (create_vector_functions, shift_fields, velocity_change,
                            get_local, set_local, dolfin_log_level)
from .field_inspector import FieldInspector
from .ufl_transformers import is_zero_ufl_expression, split_form_into_matrix
