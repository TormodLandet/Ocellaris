import re
from . import report_error

# Some imports that are useful in the code to be run
# Note the order. Dolfin overwrites NumPy which overwrites 
# the standard library math module
import math, numpy, dolfin
from math import *
from numpy import *
from dolfin import *

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
        multiline = self._validate_code(code_string)
        
        filename = '<input-file-code %s>' % description
        self.code = compile(code_string, filename, 'exec' if multiline else 'eval')
        self.is_multiline = multiline
    
    def _validate_code(self, code_string):
        """
        Check that the code is either a single expression or a valid
        multiline expression that defines the variable varname 
        """
        lines = code_string.split('\n')
        multiline = len(lines) > 1
        
        if multiline and self.var_name is not None:
            vardef = r'.*(^|\s)%s\s*=' % self.var_name
            if re.match(vardef, code_string) is None:
                report_error('Invalid: %s' % self.description,
                             'Multi line expression must define the variable "%s"'
                             % self.var_name)
        
        return multiline
    
    def run(self, **kwargs):
        """
        Run the code
        """
        # Make sure the simulation data is available 
        locals().update(self.simulation.data)
        simulation = self.simulation
        t = time = simulation.time
        it = timestep = simulation.timestep
        
        # Make sure the keyword arguments accessible
        locals().update(kwargs)
        
        if self.is_multiline:
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
