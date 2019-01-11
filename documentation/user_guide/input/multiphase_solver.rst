.. _inp_multiphase_solver:

Multi phase solver
==================

If you are creating a two fluid simulation you will have to specify some
parameters of the multi-phase solver. The best tested options are the
``SinglePhase`` and the ``BlendedAlgebraicVOF`` solvers. You can also add
your own multi phase model, see :ref:`custom_multiphase_model`.


Single phase
------------

The default value of the multi phase solver setup is to use a single phase
method:

.. code-block:: yaml

    multiphase_solver:
        type: SinglePhase

Since this is the default, the ``multiphase_solver`` section can be omitted
from the input file when running single phase simulations.


Algebraic Volume of Fluid method
--------------------------------

This is an algebraic VOF method where the interface is captured as the 0.5
level set of a colour function which should be strictly within the range
from 0 to 1 inclusive.

When using the multi phase VOF solver by specifying
:code:`type: BlendedAlgebraicVOF` the following parameters can be specified:

.. describe:: function_space_colour

    CG for continuous Galerkin, DG for discontinuous Galerkin (default)

.. describe:: polynomial_degree_colour

    The degree of the approximating polynomials, default 0 (piecewise constant)

.. describe:: num_subcycles

    Number of times the VOF method is run per time step. Running VOF with a
    smaller time step can be beneficial to decrease the Courant number. The
    VOF calculations are typically much faster than the Navier-Stokes solver.
    Default value 1, using about 5 sub cycles can to give good results in
    many cases without much impact on running time.

.. describe:: force_static

    Do not move the free surface. Can be useful for testing in some cases.

.. describe:: project_uconv_dgt0

    Project the velocity that will convect the colour function into a DGT0
    field, piecewise constant on each facet. This is consistent with how the
    finite volume method handles advection, it can be more stable (lower local
    Courant numbers) and it is mass conserving in DG0 space. Default on.

In addition you will have to specify a convection scheme for the VOF colour
function in order to keep the free surface sharp. For specifying the convection
scheme, see :ref:`inp_convection`.


Other methods
-------------

There are some other options that exist in the code base, but these are not
well tested and very likely need some work before they can run. They are

* VariableDensity (for mixtures without sharp interfaces)

  This model was actually used a bit some time ago, so it may work with a few
  modifications. It can be used with higher order descriptions of the density,
  and with an appropriate slope limiter it could potentially deal with sharp
  interface problems as well as diffusive mixes. There is a convergence test
  in ``cases/convergece-variable-density-disk`` that is not tested by default
  currently, but has been used to show that the multi phase model setup can
  produce the correct order of convergence for velocity, pressure and density.
  If you want to resurrect this method then starting with getting the
  convergence case running is probably the best option. It probably works with
  few or no code modifications.

* HeightFunction

  Implemented for 2D simuations. A single valued height function separates two
  fluid domains. Uses the height function to compute a VOF field which reuses
  most of the ``BlendedAlgebraicVOF`` machinery. Not used in a long while, so
  the code is probably somewhat bitrotted (probably no parallel support etc).

* HeightFunctionALE

  Implemented for 2D simuations. The mesh is moved according to the vertical
  fluid velocity at the free surface after each time step. Not used in a long
  while, probably somewhat bitrotted. Would need stablisation/smoothing to
  combat sawtooth instabilities if it was ever to be used for something.

  Not used in a long while, so the code is certainly bitrotted and the ALE code
  in Ocellaris in general is not tested at all at the moment, so expect
  problems. There is also probably no support for running in parallel with MPI
  etc.

* Lagrangian

  A purely Lagrangian multiphase model. The mesh is moved according to the
  calculated fluid velocity after each time step. This will obviously distort
  the mesh in allmost all calculations.

  This was implemented as a stepping stone to ALE, and to test hydrostatic
  pressure calculations where the correct answer is zero velocity everywhere
  for all time and ALE should not be necessary.

  To initialise the multi phase field the colour function must be specified in
  the input file (as initial condition for "cp"). The colour function is unity
  when rho=rho0 and nu=nu0 and zero when rho=rho1 and nu=nu1
