.. _inp_time:

Timestepping
============

This input file section sets the start time, end time, and time step. The time
step can be altered in user coding, see below. Changing the time step increment
will cause the next time step to use a first order time integration method
(Euler) to avoid consistency errors in solvers which support this (e.g. the
IPCS-A solver).

.. code-block:: yaml

    time:
        dt: 0.01
        tmax: 60.0

.. describe:: dt

    Time step

.. describe:: tstart

    Starting time, default 0.0

.. describe:: tmax

    Maximum simulation time, simulation will stop when t > tmax


Adaptive time stepping
----------------------

Adaptive time stepping can be implemented by use of hooks in the input file,
see :ref:`inp_hooks` for details.

A very simple example:

.. code-block:: yaml

    hooks:
        pre_timestep:
        -   name: decrease time step
            code: |
                if t > 10:
                    simulation.input.set_value('time/dt', 0.005)

An example of adapting the time step depending on the Courant number:

.. code-block:: yaml

    hooks:
        post_timestep:
        -   name: Adaptive time step
            code: |
                # Configuration
                TARGET_CO = 0.3
                MINIMUM_DT = 0.0001

                # Get the cell based Courant number for all previous time steps
                all_Co = simulation.reporting.timestep_xy_reports.get('Co', [0])

                # Courant number for the last time step
                Co = all_Co[-1]

                # Maximum Courant number in the last 10 time steps
                Co_window = numpy.max(all_Co[-10:])

                # Compute new dt based on the Courant numbers
                new_dt = None
                if Co > TARGET_CO:
                    new_dt = dt / 2
                elif Co_window < TARGET_CO / 6:
                    new_dt = dt * 2

                # Change the dt if required (and the new dt is not below minimum)
                if new_dt is not None and new_dt >= MINIMUM_DT:
                    simulation.log.info('Changing dt from %.5f to %.5f' % (dt, new_dt))
                    simulation.input.set_value('time/dt', new_dt)
