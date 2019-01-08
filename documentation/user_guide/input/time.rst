.. _inp_time:

Timestepping
============

This section sets the end time and time step. Currently only fixed time step is
available, though the time step can be altered in user coding. Changing the
time step increment will cause the next time step to use a first order to avoid
consistency errors in solvers with supports this (this includes the IPCS-A
solver).

.. code-block:: yaml

    time:
        dt: 0.01
        tmax: 60.0

Example user code that changes the time step. See details under
:ref:`inp_hooks`.

.. code-block:: yaml

    hooks:
        pre_timestep:
        -   name: decrease time step
            code: |
                if t > 10:
                    simulation.input['time']['dt'] = 0.005
