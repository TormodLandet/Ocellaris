.. _inp_physical_properties:

Physical properties
===================

Specify physical constants. How this is done depends mostly on your choice of
multi phase solver (see :ref:`inp_multiphase_solver`). A single phase example:

.. code-block:: yaml

    physical_properties:
        g: [0, 0, 0]
        nu: 0.001
        rho: 1.0

.. describe:: g

    The acceleration of gravity given as a list of numbers. The length of the
    list must match the number of spatial directions, e.g. 2 or 3.
    Use ``[0, -9.81]`` in 2D and ``[0, 0, -9.81]`` in 3D for "standard" gravity.


Single phase properties
-----------------------

.. describe:: nu

    The kinematic viscosity

.. describe:: rho

    The density of the fluid, defaults to ``1.0``.


VOF two phase properties
------------------------

.. describe:: nu0, rho0

    The kinematic viscosity and density of fluid 0

.. describe:: nu1, rho1

    The kinematic viscosity and density of fluid 1

For a water/air simulation fluid 0 is typically water and corresponds to
VOF colour function value ``1.0`` while fluid 1 is typically air and
corresponds to VOF colour function value ``0.0``.


Variable density properties
---------------------------

.. describe:: nu

    The kinematic viscosity of both fluids (single value)

.. describe::  rho_min, rho_max

    The range of allowable densities. Give one number for each of these settings.

