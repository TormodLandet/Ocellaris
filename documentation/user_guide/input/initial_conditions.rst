.. _inp_initial_conditions:

Initial conditions
==================

In the lid driven cavity test case both the velocity and the pressure fields
start from zero, so no initial values need to be given. The following is an
example of how to specify initial values for the Taylor-Green vortex on a 2D
square with side lengths equal to 2.0:

.. code-block:: yaml

    initial_conditions:
        up0:
            cpp_code: -sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)
        up1:
            cpp_code:  sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)
        p:
            cpp_code: -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*t)/4

.. csv-table::
   :header: "key", "Default value", "Description"

    "initial_conditions/var_name/cpp_code", "**required input**", "C++ code that gives the value of the field at each point. Variables ``rho``, ``nu`` and ``t`` are available"
