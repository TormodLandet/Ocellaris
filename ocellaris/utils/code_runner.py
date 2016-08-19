import re
from . import ocellaris_error, run_debug_console

# Some imports that are useful in the code to be run
# Note the order. Dolfin overwrites NumPy which overwrites 
# the standard library math module
import math, numpy, dolfin
from math import *
from numpy import *
from dolfin import *


__all__ = ['RunnablePythonString', 'CodedExpression']


class RunnablePythonString(object):
    def __init__(self, simulation, code_string, description, var_name=None):
        """
        This class handles Python code that is given on the
        input file
        
        It the code contains a newline it must be the core
        of a function and is run with exec() otherwise it is
        assumed to be an expression and is run with eval()
        
        If varname is specified then any multiline code block
        must define this variable
        """
        self.simulation = simulation
        self.description = description
        self.var_name = var_name
        needs_exec = self._validate_code(code_string)
        
        filename = '<input-file-code %s>' % description
        self.code = compile(code_string, filename, 'exec' if needs_exec else 'eval')
        self.needs_exec = needs_exec
    
    def _validate_code(self, code_string):
        """
        Check that the code is either a single expression or a valid
        multiline expression that defines the variable varname 
        """
        # Does this code define the variable, var_name = ...
        # or assign to an element, var_name[i] = ... ?
        if self.var_name is not None:
            vardef = r'.*(^|\s)%s(\[\w\])?\s*=' % self.var_name
            has_vardef = re.search(vardef, code_string) is not None
        else:
            has_vardef = False
        
        multiline = '\n' in code_string
        needs_exec = multiline or has_vardef
        
        if needs_exec and self.var_name is not None and not has_vardef:
            ocellaris_error('Invalid: %s' % self.description,
                            'Multi line expression must define the variable "%s"'
                            % self.var_name)
        
        return needs_exec
    
    def run(self, **kwargs):
        """
        Run the code
        """
        # Make sure the simulation data is available
        simulation = self.simulation 
        locals().update(simulation.data)
        t = time = simulation.time
        it = timestep = simulation.timestep
        dt = simulation.dt
        ndim = simulation.ndim
        
        # Make sure the keyword arguments accessible
        locals().update(kwargs)
        
        # Make sure the user constants are accessible
        user_constants = simulation.input.get_value('user_code/constants', {}, 'dict(string:float)')
        constants = {}
        for name, value in user_constants.iteritems():
            constants[name] = value
        locals().update(constants)
        
        if self.needs_exec:
            exec(self.code)
            if self.var_name is not None:
                # The code defined a variable. Return it
                return locals()[self.var_name]
            else:
                # No return value
                return
        else:
            # Return the result of evaluating the expression
            return eval(self.code)


def CodedExpression(simulation, code_string, description, value_shape=()):
    """
    This Expression sub-class factory creates objects that run the given
    RunnablePythonString object when asked to evaluate
    """
    # I guess dolfin overloads __new__ in some strange way ???
    # This type of thing should really be unnecessary ...
    if value_shape == ():
        expr = CodedExpression0()
    elif value_shape == (2,):
        expr = CodedExpression2()
    elif value_shape == (3,):
        expr = CodedExpression2()
    
    expr.runnable = RunnablePythonString(simulation, code_string, description, 'value')
    return expr


################################################################################
# We need to subclass once per value_shape() for some reason

class CodedExpression0(dolfin.Expression):
    def eval_cell(self, value, x, ufc_cell):
        self.runnable.run(value=value, x=x, ufc_cell=ufc_cell)


class CodedExpression2(dolfin.Expression):
    def eval_cell(self, value, x, ufc_cell):
        self.runnable.run(value=value, x=x, ufc_cell=ufc_cell)
    
    def value_shape(self):
        return (2,)


class CodedExpression3(dolfin.Expression):
    def eval_cell(self, value, x, ufc_cell):
        self.runnable.run(value=value, x=x, ufc_cell=ufc_cell)

    def value_shape(self):
        return (3,)
