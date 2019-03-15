# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import re
from . import ocellaris_error, run_debug_console

# Some imports that are useful in the code to be run
# Note the order. Dolfin overwrites NumPy which overwrites
# the standard library math module
import math
import numpy
import dolfin
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
            ocellaris_error(
                'Invalid: %s' % self.description,
                'Multi line expression must define the variable "%s"' % self.var_name,
            )

        return needs_exec

    def run(self, **kwargs):
        """
        Run the code
        """
        # Make sure some useful variables are available
        sim = simulation = self.simulation
        t = time = simulation.time
        it = timestep = simulation.timestep
        dt = simulation.dt
        ndim = simulation.ndim

        # Make the simulation data is available
        code_locals = locals()
        code_locals.update(simulation.data)

        # Make sure the keyword arguments accessible
        code_locals.update(kwargs)

        # Make sure the user constants are accessible
        user_constants = simulation.input.get_value(
            'user_code/constants', {}, 'dict(string:basic)'
        )
        constants = {}
        for name, value in user_constants.items():
            constants[name] = value
        code_locals.update(constants)

        if self.needs_exec:
            try:
                exec(self.code, globals(), code_locals)
            except Exception:
                sim.log.error('Python code exec failed for the below code:')
                sim.log.error(self.code)
                raise

            if self.var_name is not None:
                # The code defined a variable. Return it
                return code_locals[self.var_name]
            else:
                # No return value
                return
        else:
            # Return the result of evaluating the expression
            try:
                return eval(self.code, globals(), code_locals)
            except Exception:
                sim.log.error('Python code eval failed for the below code:')
                sim.log.error(self.code)
                raise


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

    expr._runnable = RunnablePythonString(simulation, code_string, description, 'value')
    return expr


################################################################################
# We need to subclass once per value_shape() for some reason


class CodedExpression0(dolfin.UserExpression):
    def eval_cell(self, value, x, ufc_cell):
        ret = self._runnable.run(value=value, x=x, ufc_cell=ufc_cell)
        if not self._runnable.needs_exec:
            value[:] = ret


class CodedExpression2(dolfin.UserExpression):
    def eval_cell(self, value, x, ufc_cell):
        ret = self._runnable.run(value=value, x=x, ufc_cell=ufc_cell)
        if not self._runnable.needs_exec:
            value[:] = ret

    def value_shape(self):
        return (2,)


class CodedExpression3(dolfin.UserExpression):
    def eval_cell(self, value, x, ufc_cell):
        ret = self._runnable.run(value=value, x=x, ufc_cell=ufc_cell)
        if not self._runnable.needs_exec:
            value[:] = ret

    def value_shape(self):
        return (3,)
