Advanced topics
===============

You can extend Ocellaris in more ways than just adding user code for
customizing the existing fluid solvers. You can for instance add brand new
solvers while still taking advantage of the Ocellaris framework for input and
output handling, post-processing and more.

.. contents:: Contents
    :local:


.. _custom_solver:

Writing a custom solver
-----------------------

The architecture of Ocellaris is quite customizable, so if you want to add a
custom equation solver, multiphase model, convection scheme, slope limiter or
more it is possible to do so without touching the Ocellaris code, but by
crafting a special input file.

A minimal custom solver must implement the following:

.. code-block:: python

    import dolfin
    from ocellaris.solvers import Solver, register_solver

    @register_solver('MyCustomSolver')
    class ACustomSolverClass(Solver):
        description = "This is my custom solver!"

        @classmethod
        def create_function_spaces(cls, simulation):
            mesh = simulation.data['mesh']
            simulation.data['Vmyfunc'] = dolfin.FunctionSpace(mesh, family, degree)

        def __init__(self, simulation):
            self.simulation = simulation

        def run(self):
            sim = self.simulation
            sim.hooks.simulation_started()

            dt = 0.1
            t = 0
            for i in range(1, 10):
                t += dt
                sim.hooks.new_timestep(timestep_number=i, t=t, dt=dt)

                #  .... do something ...

                sim.hooks.end_timestep()
            sim.hooks.simulation_ended(success=True)

Ocellaris will first run ``create_function_spaces`` and then a bit later the
``__init__`` method will run. A full example is given next, a custom solver
for solving a Poisson equation for an unknown function, phi, using the
discontinuous Galerkin symmetric interior penalty method. This file can be
found in ``custom_solver.py`` under the ``demos`` folder in the Ocellaris
source code.

.. literalinclude:: ../../demos/custom_solver.py
   :language: python

The input file that goes along with this solver makes sure that the solver
Python file is loaded and the sets up the boundary conditions for the
``phi`` field.

.. literalinclude:: ../../demos/custom_solver.inp
   :language: yaml


.. _custom_boundary_condition:

Writing a custom boundary condition
-----------------------------------

Here the procedure to create a custom Dirichlet boundary condition is sketched.
Creating custom Neuman or Robin BCs (Boundary Conditions) follows roughly the
same setup, look at the implementation of those (or FreeSlip etc) in the source
code for pointers to what steps must be taken.

The Ocellaris input file must load the custom BC Python file to activate it,
and then the input file can reference the custom BC type later on:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0

    user_code:
        modules:
        -   custom_bc_module

    # ...

    boundary_conditions:
    -   name: all walls
        selector: code
        inside_code: on_boundary
        u0:
            type: MyCustomBC

The ``custom_bc_module.py`` file defines ``MyCustomBC``. The code can be
something like this (untested, but the overall design should work):

.. code-block:: python

    import dolfin
    from ocellaris.solver_parts.boundary_conditions.dirichlet import (
        register_boundary_condition,
        BoundaryConditionCreator,
        OcellarisDirichletBC)

    @register_boundary_condition('MyCustomBC')
    class CustomDirichletBoundary(BoundaryConditionCreator):
        description = 'A custom Dirichlet boundary condition'

        def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
            """
            Dirichlet boundary condition that does X
            """
            self.simulation = simulation

            # This specific BC only works for the velocity component "u0"
            # (this is not very realistic, but it simplifies this example)
            assert var_name == 'u0'
            V = simulation.data['Vu']

            # Define a function that holds the BC value
            bc_val = dolfin.Function(V)

            bc = OcellarisDirichletBC(simulation, V, bc_val, subdomains, subdomain_id)
            bcs = simulation.data['dirichlet_bcs']
            bcs.setdefault(var_name, []).append(bc)
            self.simulation.log.info('    Custom BC for %s' % var_name)

            # Update this BC before each time step
            self.bc_val_func = bc_val
            simulation.hooks.add_pre_timestep_hook(self.update, 'CustomBC update')

        def update(self, timestep_number, t, dt):
            """
            This code runs before the Navier-Stoke solver each time step and
            can change the BC value that is used in the assembly of the system
            matrices. It must modify the function bc_val that was sent to the
            OcellarisDirichletBC class above.
            """
            some_code_to_update_this_bc(self.bc_val_func)


.. _custom_multiphase_model:

Writing a custom multiphase model
---------------------------------

Connecting a custom multi phase model to Ocellaris is done in the same way as
connecting the main equation solver described above. How to register the model
is shown in the following snippet

.. code-block:: python

    from ocellaris.solver_parts.multiphase import MultiPhaseModel, \
                                                  register_multi_phase_model

    class MyCustomMultiPhaseModel(MultiPhaseModel):
        description = 'A fantastic multiphase model!'

        # See the SinglePhase model for an example of the required
        # API interface needed to work inside Ocellaris


Writing other custom parts
--------------------------

Many more items are pluggable, some examples are

* Convection schemes (used mostly in the VOF multi phase solver)
* Known fields (for initial and boundary conditions, forcing zones and more)
* Slope limiters (scalar fields)
* Velocity slope limiters (solenoidal vector fields)

You will have to look in the source code for the required API to implement, but
the method for connecting the custom code to Ocellaris is the same as what is
described above for the main equation solver. The registering functions are:

.. code-block:: python

    from ocellaris.solver_parts.convection import register_convection_scheme
    from ocellaris.solver_parts.fields import register_known_field
    from ocellaris.solver_parts.slope_limiter import register_slope_limiter
    from ocellaris.solver_parts.slope_limiter_velocity import register_velocity_slope_limiter

If you want to create a custom forcing zone or something else that is not
currently pluggable, then please create a feature request on the issue tracker
or submit a pull request with code following the same general framework as the
pluggable pieces above.
