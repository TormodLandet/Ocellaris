import dolfin
from . import report_error

def make_expression(simulation, cpp_code, description):
    """
    Create a C++ expression with parameters like time and all scalars in 
    simulation.data available (nu and rho for single phase simulations) 
    """
    available_vars = get_vars(simulation)    
    try:
        return dolfin.Expression(cpp_code, **available_vars)
    except Exception as e:
        vardesc = '\n  - '.join('%s (%s)' % (name, type(value)) for name, value in available_vars.items())
        errormsg  = str(e)
        report_error('Error in C++ code',
                     'The C++ code for %s does not compile.'
                     '\n\nCode:\n%s'
                     '\n\nGiven variables:\n  - %s'
                     '\n\nError:\n%s' % (description, cpp_code, vardesc, errormsg))


def get_vars(simulation):
    """
    Make a dictionary of variables to send to the C++ expression. Returns the
    time "t" and any scalar quantity in simulation.data
    """
    available_vars = {'t': simulation.time, 'it': simulation.timestep}
    for name, value in simulation.data.items():
        if isinstance(value, (float, int, long)):
            available_vars[name] = value
        elif isinstance(value, dolfin.Constant) and value.ufl_shape == ():
            available_vars[name] = value
    
    # Sanity check of variable names
    for name in available_vars:
        assert not hasattr(dolfin.Expression, name)
    
    return available_vars


def ocellaris_project(simulation, cpp_code, description, V, function=None):
    """
    Create a C++ expression with parameters like time and all scalars in 
    simulation.data available (nu and rho for single phase simulations) 
    
    Project the expression into a dolfin.Function. The function is either
    provided or a function will be created 
    """
    # Compile the C++ code
    expr = make_expression(simulation, cpp_code, description)
    
    # Project to the function space and return the projected function
    if function is None:
        function = dolfin.Function(V)
    dolfin.project(expr, V=V, function=function)
    return function


def OcellarisCppExpression(simulation, cpp_code, description, update=False):
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
        # Update the expression with new values of time an similar
        # scalar quantities
        available_vars = get_vars(simulation)
        for name, value in available_vars.items():
            if hasattr(expression, name):
                setattr(expression, name, value)
    
    # Create the expression
    expression = make_expression(simulation, cpp_code, description)
    
    # Return the expression. Optionally register an update each time step
    if update:
        simulation.hooks.add_pre_timestep_hook(updater, 'Update C++ expression "%s"' % description)
    
    return expression
