# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import traceback
import dolfin
from . import ocellaris_error, timeit


def make_expression(simulation, cpp_code, description, element, params=None):
    """
    Create a C++ expression with parameters like time and all scalars in
    simulation.data available (nu and rho for single phase simulations)
    """
    if isinstance(cpp_code, (float, int)):
        cpp_code = repr(cpp_code)

    if isinstance(element, int):
        degree = element
        element = None
    else:
        degree = None

    # Get variables that can be used in the expression, but are defined
    # elsewhere, either in Ocellaris (timestep, density etc), or defined
    # in the user constants section in the input file.
    available_vars = get_vars(simulation)
    if params:
        available_vars.update(params)

    # Get all identifiers found in the C++ code. Variables from the
    # available_vars dictionary that are not found as identifiers in
    # the code are not needed and can be dropped. This simplifies and
    # potentially speeds up the generated code
    identifiers = get_identifiers(cpp_code)
    available_vars = {k: v for k, v in available_vars.items() if k in identifiers}

    try:
        return dolfin.Expression(
            cpp_code, element=element, degree=degree, **available_vars
        )
    except Exception:
        vardesc = '\n  - '.join(
            '%s (%s)' % (name, type(value)) for name, value in available_vars.items()
        )
        errormsg = traceback.format_exc()
        ocellaris_error(
            'Error in C++ code',
            'The C++ code for %s does not compile.'
            '\n\nCode:\n%s'
            '\n\nGiven variables:\n  - %s'
            '\n\nError:\n%s' % (description, cpp_code, vardesc, errormsg),
        )


def get_vars(simulation):
    """
    Make a dictionary of variables to send to the C++ expression. Returns the
    time "t" and any scalar quantity in simulation.data
    """
    available_vars = {
        't': simulation.time,
        'dt': dolfin.Constant(simulation.dt),
        'it': simulation.timestep,
        'ndim': simulation.ndim,
    }

    # Simulation fields etc
    for name, value in simulation.data.items():
        if isinstance(value, (float, int)):
            available_vars[name] = value
        elif isinstance(value, dolfin.Constant) and value.ufl_shape == ():
            available_vars[name] = value

    # User constants
    user_constants = simulation.input.get_value(
        'user_code/constants', {}, 'dict(string:basic)'
    )
    for name, value in user_constants.items():
        if isinstance(value, (int, float)):
            available_vars[name] = value

    # Sanity check of variable names
    for name in available_vars:
        assert not hasattr(dolfin.Expression, name)

    return available_vars


def get_identifiers(code_string):
    """
    Return all valid identifiers found in the C++ code string by finding
    uninterrupted strings that start with a character or an underscore and
    continues with characters, underscores or digits. Any spaces, tabs,
    newlines, paranthesis, punctuations or newlines will mark the end of
    an identifier (anything that is not underscore or alphanumeric).

    The returned set will include reserved keywords (if, for, while ...),
    names of types (float, double, int ...), names of functions that are
    called (cos, sin ...) in addition to any variables that are used or
    defined in the code.
    """
    identifiers = set()
    identifier = []
    for c in code_string + ' ':
        if identifier and (c == '_' or c.isalnum()):
            identifier.append(c)
        elif c == '_' or c.isalpha():
            identifier = [c]
        elif identifier:
            identifiers.add(''.join(identifier))
            identifier = None
    return identifiers


@timeit
def ocellaris_interpolate(
    simulation, cpp_code, description, V, function=None, params=None
):
    """
    Create a C++ expression with parameters like time and all scalars in
    simulation.data available (nu and rho for single phase simulations)

    Interpolate the expression into a dolfin.Function. The results can be
    returned in a provided function, or a new function will be returned
    """
    # Compile the C++ code
    with dolfin.Timer('Ocellaris make expression'):
        expr = make_expression(
            simulation, cpp_code, description, element=V.ufl_element(), params=params
        )

    if function is None:
        function = dolfin.Function(V)
    else:
        Vf = function.function_space()
        if (
            not Vf.ufl_element().family() == V.ufl_element().family()
            and Vf.dim() == V.dim()
        ):
            ocellaris_error(
                'Error in ocellaris_interpolate',
                'Provided function is not in the specified function space V',
            )

    with dolfin.Timer('Ocellaris interpolate expression'):
        return function.interpolate(expr)


def OcellarisCppExpression(
    simulation,
    cpp_code,
    description,
    degree,
    update=True,
    return_updater=False,
    params=None,
    quad_degree=None,
):
    """
    Create a dolfin.Expression and make sure it has variables like time
    available when executing.

    If update is True: all variables are updated at the start of each time
    step. This is useful for boundary conditions that depend on time
    """

    def updater(timestep_number, t, dt):
        """
        Called by simulation.hooks on the start of each time step
        """
        # Update the expression with new values of time and similar
        # scalar quantities
        available_vars = get_vars(simulation)
        for name, value in available_vars.items():
            if name in expression.user_parameters:
                expression.user_parameters[name] = value

    if quad_degree is not None:
        mesh = simulation.data['mesh']
        element = dolfin.FiniteElement(
            'Quadrature', mesh.ufl_cell(), quad_degree, quad_scheme='default'
        )
    else:
        # Element = integer creates a CG element with the given degree
        element = degree

    # Create the expression
    expression = make_expression(simulation, cpp_code, description, element, params)

    # Return the expression. Optionally register an update each time step
    if update:
        simulation.hooks.add_pre_timestep_hook(
            updater, 'Update C++ expression "%s"' % description, 'Update C++ expression'
        )

    if return_updater:
        return expression, updater
    else:
        return expression
