.. _inp_solver:

Navier Stokes solver
====================

Some parameters are shared between all the available velocity-pressure solvers.
All the parameters in the following example have sensible defaults except for
the solver ``type`` which you must set. For the other parameters in the example
the values shown are the default values. The possible values for solver type
(IPCS-A, SIMPLE, PISO etc) are listed in the sections below with a brief
description.

.. code-block:: yaml

    solver:
        type: IPCS-A
        num_inner_iter: 10
        allowable_error_inner: 1.0e-10
        polynomial_degree_pressure: 1
        polynomial_degree_velocity: 2
        function_space_pressure: DG
        function_space_velocity: DG
        u:
            # see linear solver documentation below
        p:
            # see linear solver documentation below

The inner iterations (pressure correction iterations) will run a maximum of
``num_inner_iter`` times for each time step, but the iterations will exit early
if the :math:`l^2` error of the difference between the predicted and corrected
velocity field is less than the given value ``allowable_error_inner``.

Some control parameters exist outside the common ones shown above, but none of
these are of the type that a normal user would probably need to change, so they
are only documented in the source code of the individual solvers.

The following parameters are relevant for under-relaxed solver implementations
(SIMPLE, PISO, PIMPLE):

.. describe::  relaxation_u, relaxation_p

    Relaxation factors. A value of 1.0 means no relaxation, 0.0 means no update
    at all (pointless). A value of 0.5 means that the result is an even blend
    of the computed value and the previous iteration value

.. describe::  relaxation_u_last_iter, relaxation_p_last_iter

    Some solvers will differentiate the last inner iteration from all other
    iterations. These parameters default to 1.0 in order to perform a "propper"
    update at the end of a time step with no relaxation applied.


IPCS-A
------

Incremental Pressure Correction Scheme on Algebraic form. This is an iterative
Chorin/Temam type pressure correction solver.


IPCS-D
------

Incremental Pressure Correction Scheme on Differential form. This is an
iterative Chorin/Temam type pressure correction solver where the pressure
correction Poisson equation is assembled from an elliptic operator and not
algebraicly from matrices. The divergence of the velocity field is hence not
very low and the method is not so strongly recommended for DG FEM, but it is
one of the most common solvers for the Navier-Stokes equations outside of DG
FEM and it has a smaller numerical stencil and may be faster than the IPCS-A
method.


SIMPLE
------

Semi-Implicit Method for Pressure-Linked Equations. The implementation of the
algorithm is based on Klein, Kummer, Keil & Oberlack (2015).


PISO
----

The pressure correction method by Issa (1986), Pressure-Implicit with Splitting
of Operators. PISO adds an additional correction step to the SIMPLE algorithm.


PIMPLE
------

A Navier-Stokes solver based on the PIMPLE algorithm as implemented in OpenFOAM
and partially described in the PhD thesis of Jasak (1996; the PISO loop only).

.. describe::  num_pressure_corr

    The number of PISO iterations for each PIMPLE loop (the number of PIMPLE
    loops is controlled by the standard ``num_inner_iter`` parameter).


Coupled
-------

Solves the velocity-pressure saddle point block-matrix equation system coupled.
Do not use this solver for large meshes. Even when using the multi-cpu
distributed multi frontal MUMPS or SuperLU_dist direct solvers there is a quite
small (perhaps around 1 million on a recent workstation?) limit to how many
degrees of freedom can be computed. For very small examples it may be faster
than using pressure-correction iterations and there is no resulting splitting
error which makes it great for testing and benchmarking the split solvers.

No block-system preconditioners are available in Ocellaris for the coupled
Navier-Stokes solver, so iterative linear solvers will either not converge or
perhaps "converge" to nonsensical solutions. Only use with direct solvers!


Analytical
----------

Use the initial condition C++ code (possibly containing the time variable ``t``
which will be updated for each time step) to define the velocity and pressure
for all time steps. This can be usefull for testing other parts of the
Ocellaris solution framework with a known Navier-Stokes solution.


Specifying the linear solver
----------------------------

All equation systems that require global solves, like the velocity, pressure
and potentially multi phase models, will have their own optional definition of
the linear solver. These can be described in two ways, the simple FEniCS DOLFIN
based setup where some limited configuration is possible, or the full PETSc KSP
setup where all of the PETSc options are configurable plus a few options added
by Ocellaris.

It is recommended to use the KSP setup. It is the default, it is more powerfull
and it can do everything supported by the FEniCS DOLFIN setup. The DOLFIN setup
is kept for comparison and to be able to test the exact same setup used by
"normal" FEniCS codes.


PETSC KSP solver setup (use_ksp = yes)
......................................

This linear solver setup is used by most linear solvers inside Ocellaris. Most
solvers set reasonable defaults. Use these as starting points for your own
experimentations. The Ocellaris log file shows the setup which is used for the
different linear solvers in your simulation.

.. code-block:: yaml

    solver:
        u:
            use_ksp: yes
            petsc_ksp_type: gmres
            petsc_pc_type: asm
            petsc_ksp_initial_guess_nonzero: yes
            inner_iter_rtol: [1.0e-15, 1.0e-15, 1.0e-15]
            inner_iter_atol: [1.0e-15, 1.0e-15, 1.0e-15]
            inner_iter_max_it: [100, 100, 100]

.. describe:: use_ksp: yes

    Signal that we want to use the KSP solver setup (this is default in most
    situations).

.. describe:: petsc_XXXX

    Any PETSc parameter. Examples: ``ksp_type`` sets the solver name and
    ``pc_type`` sets the preconditioner name. Look at the PETSc documentation
    for the full list of tunable parameters, or give ``petsc_help: 'ENABLED'``
    to get a dump of possible parameters (the program will exit after giving
    the parameter listing).

.. describe:: inner_iter_control

    The number of iterations and tolerances in the Krylov solver can be set for
    three categories of solves. The first X inner iterations (pressure
    correction iterations in the Navier-Stokes solver), the last Y inner
    iterations and the rest of the iterations (the middle number). The numbers
    X and Y are set by ``inner_iter_control: [X, Y]``. The default values are
    ``X=Y=3``.

.. describe:: inner_iter_rtol, inner_iter_atol, inner_iter_max_it

    The relative and absolute tolerances in the Krylov solver (default values
    are typically ``rtol = 1.0e-10`` and ``atol = 1.0e-15``). The maximum
    number of Krylov iterations is by default ``100`` for most solvers. If the
    solution is not converged the procedure will just continue, it is not
    always necessary to fully converge when applying an iterative solver, at
    least not in the inner first iterations (see below note on iterations).

.. note::

    Inner iterations refer to the main iterations inside each time step,
    typically pressure correction iterations (implemented in code inside
    Ocellaris). Krylov iterations refer to iterations inside the linear
    equation solver (provided by PETSc). The Krylov iterations are nested
    inside the inner iterations which are nested inside the time loop.


FEniCS DOLFIN solver setup (use_ksp = no)
.........................................

.. code-block:: yaml

    solver:
        u:
            use_ksp: no
            solver: gmres
            preconditioner: additive_schwarz
            parameters:
                any_parameter_supported_by_dolfin: valid_value

.. describe:: use_ksp: no

    Signal that we want to use the simplified setup

.. describe::  solver, preconditioner

    The names of the preconditioner and linear solver. Any values (string)
    supported by FEniCS DOLFIN are supported. The default values in FEniCS
    are used if none are specified (bad idea for large systems)

.. describe::  parameters

    Any parameter keys and values supported by FEniCS DOLFIN. See the DOLFIN
    documentation for these.
