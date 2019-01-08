.. _inp_user_code:

User constants and code
=======================

You can specify constants that can be used in subsequent sections to make
the input file easily configurable. You can also specify some code that
will run right after the input file has been read, before any of the 
simulation setup such as loading the mesh has been done. You can even
change the input by accessing the ``simulation.input`` object since no
parts of Ocellaris has accessed the input yet.

.. code-block:: yaml

    user_code:
        constants:
            L: 200       # channel length
            theta: 30    # angle
        code: |
            import subprocess
            subprocess.call(['command', 'to', 'generate', 'mesh'])

Example of using the constants in later sections of the input file:

.. code-block:: yaml

    some:
        section:
            param1: 4.3
            param2: py$ 2.3 * L * sin(theta)
            cpp_code: 'x[0] + L * sin(theta)' 

Any value (except inside the ``user_code/constants`` block) can be given as
a string starting with ``py$``. Ocellaris will then execute the given Python
code to produce the value to be used in Ocellaris just as if you had written
the value directly into the input file. The Python code you give can evaluate
to a list, string, number...

Code given as strings in the input file, either Python or C++ can also use
the constants as can be seen in the example. These are typically expressions
defining initial or boundary values. You can even combine these functions:

.. code-block:: yaml

    some-section:
        cpp_code: py$ 'x[0] + L * sin(theta)'.replace('theta', 'theta + L') 

This can be handy if you give the C++ code to compute the value of a field
as a user constant string, and then you can use python code to replace the
variable  ``t`` in the string with ``(t - dt)`` in order to specify the two
initial conditions, both at ``t=0`` and ``t=0-dt`` without having to repeat
the C++ code. This can, e.g., be used to describe a Taylor-Green vortex in
such a way that the time stepping can be second order from the first time
step (normally the first time setp is first order accurate since only one
initial condition is specified:


.. code-block:: yaml

    user_code:
        constants:
            u0a: '-sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)'
            u1a: ' sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)'

    initial_conditions:
        up0:
            cpp_code: py$ u0a
        up1:
            cpp_code: py$ u1a
        upp0:
            cpp_code: py$ u0a.replace('*t)', '*(t - dt))')
        upp1:
            cpp_code: py$ u1a.replace('*t)', '*(t - dt))')
