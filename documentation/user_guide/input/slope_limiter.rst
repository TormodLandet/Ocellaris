.. inp_slope_limiter:

Slope limiter
=============

Slope limiters are needed to combat non-linear convective instabilities when
jumps or shocks are present in the solution or coefficients. Ocellaris uses a
special velocity slope limiting procedure to stabilise the momentum equation,
see :cite:`landet_slope_lim` for details.

.. describe::  method

    The method employed to perform the slope limiting

    Scalar methods:

    * None
    * OnlyBound
    * HierarchicalTaylor

    Velocity methods:

    * Componentwise
    * Solenoidal (not ported to FEniCS 2018.1)

    Default value: None

.. describe:: enforce_bounds

    Take the min and max of the initial scalar field and clamp all subsequent
    solutions to be within these global bounds. Default value: off.

.. describe:: skip_boundaries

    A list of boundary names to skip, cells sharing a facet with these
    boundaries are not limited. Default value: ``[]``.

.. describe:: plot

    Plot some information about the limiter. How much "wigglyness" is detected
    etc. Mostly for debug pruposes. Default value: off.

.. describe:: use_cpp

    Use the C++ implementation and not the Python implementation if both exist.
    Default value: on.


Options for the HierarchicalTaylor limiter
------------------------------------------

.. describe:: enforce_bcs

    Use the Dirichlet boundary conditions for the field in the limiter, not
    only to set the "safe range" of values, but also to clamp the boundary
    values directly to the BC values in the limiter. Default value: on.

.. describe:: use_weak_bcs

    Trust the solution at the Dirichlet boundaries. The weak Dirichlet BCs
    may have kept the solution at the boundary from blowing up. Use the
    solution instead of the user specified Dirichlet value. Default value: on.

.. describe:: trust_robin_dval

    The solution is expected to be close to the Robin BC dval value, so trust
    this value and use it to extend the "safe range" of values for vertices at
    the Robin boundary. Default value: on.


Options for velocity field limiters
-----------------------------------

.. describe::  comp_method

    Only relevant for ``Componentwise`` limiter. Select which scalar limiter to
    use for the velocity components. Default value: None.

.. describe::  limit_conv

    Limit the convecting velocity. Default value: off.
