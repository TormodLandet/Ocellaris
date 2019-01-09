.. _inp_initial_conditions:

Initial conditions
==================

Specify initial conditions for the simulation. You can give initial conditions
for the previous time step and the one before that if you want to start with a
higher order time stepping scheme. The following is an example of how to
specify initial values for the Taylor-Green vortex on a 2D square with side
lengths equal to 2.0:

.. code-block:: yaml

    initial_conditions:
        up0:
            cpp_code: -sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)
        up1:
            cpp_code:  sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)
        p:
            cpp_code: -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*t)/4

Previous time step values can be given for ``upX``, ``cp`` and ``p`` (in a
multi-phase VOF simulation). The pressure will only be used in pressure
correction methods to provide an initial guess. You can also give values
for the time step before the previous by specifying ``uppX`` and ``cpp``. The
"p" is short for previous (except when it means pressure).

.. describe:: cpp_code

    C++ code that gives the value of the field at each point. Variables
    ``rho``, ``nu`` and ``t`` are available in addition to the location
    ``x[i]`` and any constants given in ``user_constants``, see
    :ref:`inp_user_code`.

.. describe:: file

    You can load initial conditions from a restart file

    .. code-block:: yaml

        initial_conditions:
            file:
                h5_file: my_restart_file.h5
                same_mesh: yes

.. describe:: function

    You can give the name of a field function instead of using C++ code

    .. code-block:: yaml

        initial_conditions:
            cp:
                function: waves/c
            up0:
                function: waves/uhoriz
            up2:
                function: waves/uvert

    See :ref:`inp_fields` for more information on known fields.
