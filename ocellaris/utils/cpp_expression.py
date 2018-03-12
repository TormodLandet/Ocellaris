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
    
    available_vars = get_vars(simulation)
    if params:
        available_vars.update(params)
    
    try:
        return dolfin.Expression(cpp_code, element=element, degree=degree, **available_vars)
    except Exception:
        vardesc = '\n  - '.join('%s (%s)' % (name, type(value)) for name, value in available_vars.items())
        errormsg = traceback.format_exc()
        ocellaris_error('Error in C++ code',
                        'The C++ code for %s does not compile.'
                        '\n\nCode:\n%s'
                        '\n\nGiven variables:\n  - %s'
                        '\n\nError:\n%s' % (description, cpp_code, vardesc, errormsg))


def get_vars(simulation):
    """
    Make a dictionary of variables to send to the C++ expression. Returns the
    time "t" and any scalar quantity in simulation.data
    """
    available_vars = {'t': simulation.time,
                      'dt': dolfin.Constant(simulation.dt),
                      'it': simulation.timestep,
                      'ndim': simulation.ndim}
    
    # Simulation fields etc
    for name, value in simulation.data.items():
        if isinstance(value, (float, int)):
            available_vars[name] = value
        elif isinstance(value, dolfin.Constant) and value.ufl_shape == ():
            available_vars[name] = value
    
    # User constants
    user_constants = simulation.input.get_value('user_code/constants', {}, 'dict(string:basic)')
    for name, value in user_constants.items():
        if isinstance(value, (int, float)):
            available_vars[name] = value
    
    # Sanity check of variable names
    for name in available_vars:
        assert not hasattr(dolfin.Expression, name)
    
    return available_vars


@timeit
def ocellaris_interpolate(simulation, cpp_code, description, V,
                          function=None, params=None):
    """
    Create a C++ expression with parameters like time and all scalars in 
    simulation.data available (nu and rho for single phase simulations) 
    
    Interpolate the expression into a dolfin.Function. The results can be
    returned in a provided function, or a new function will be returned 
    """
    # Compile the C++ code
    with dolfin.Timer('Ocellaris make expression'):
        expr = make_expression(simulation, cpp_code, description,
                               element=V.ufl_element(), params=params)
    
    if function is None:
        function = dolfin.Function(V)
    else:
        Vf = function.function_space()
        if not Vf.ufl_element().family() == V.ufl_element().family() and Vf.dim() == V.dim():
            ocellaris_error('Error in ocellaris_interpolate',
                            'Provided function is not in the specified function space V')
    
    with dolfin.Timer('Ocellaris interpolate expression'):
        return function.interpolate(expr)


def OcellarisCppExpression(simulation, cpp_code, description, degree,
                           update=True, return_updater=False,
                           params=None):
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
    
    # Create the expression
    expression = make_expression(simulation, cpp_code, description, degree, params)
    
    # Return the expression. Optionally register an update each time step
    if update:
        simulation.hooks.add_pre_timestep_hook(updater, 'Update C++ expression "%s"' % description,
                                               'Update C++ expression')
    
    if return_updater:
        return expression, updater
    else:
        return expression
