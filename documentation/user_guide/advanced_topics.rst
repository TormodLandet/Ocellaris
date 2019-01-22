.. _user_guide_advanced_topics:

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


.. _custom_known_field:

Writing a custom known field
----------------------------

Known fields are Python classes that provide a ``get_function`` method. The
following example shows building and enabling a very simple field:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0

    user_code:
        modules:
        -   custom_field

    fields:
    -   name: myfield
        type: MyCustomField
        myval: 42.0


    # ...
    # input can now refer to the 'myfield/gamma' function

The ``custom_field.py`` file defines ``MyCustomField`` which defines a
``gamma`` function. The code can be something like this:

.. code-block:: python

    import dolfin
    from ocellaris.solver_parts.fields import register_known_field, KnownField

    @register_known_field('MyCustomField')
    class SimpleCustomField(KnownField):
        description = 'Simple custom field'

        def __init__(self, simulation, field_inp):
            """
            A scalar field
            """
            self.simulation = simulation
            self.func = None
            self.value = field_inp.get_value('myval', required_type='float')

        def get_variable(self, name):
            if name != 'gamma':
                ocellaris_error(
                    'Custom field does not define %r' % name,
                    'This custom field defines only gamma',
                )
            if self.func is None:
                V = self.simulation.data['Vu']
                self.func = dolfin.Function(V)
                arr = self.func.vector().get_local()
                arr[:] = self.value
                self.func.vector().set_local(arr)
                self.func.vector().apply('insert')
            return self.func


A custom wave field
...................

A custom wave field can be made in the same way as the simple field above, but
some benefits can be had from using the ``BaseWaveField`` base class which
handles partial filling of cells in VOF simulations among other things.


.. code-block:: python

    import dolfin
    from ocellaris.utils import ocellaris_error
    from ocellaris.solver_parts.fields import register_known_field
    from ocellaris.solver_parts.fields.base_wave_field import BaseWaveField


    @register_known_field('MyCustomWaves')
    class CustomWaveField(BaseWaveField):
        description = 'Custom wave field'

        def __init__(self, simulation, field_inp):
            """
            A custom wave field
            """
            super(RaschiiWaveField, self).__init__(simulation, field_inp)
            simulation.log.info('Creating a custom wave field %r' % self.name)

            # Read input
            wave_height = field_inp.get_value('wave_height', required_type='float')
            h = field_inp.get_value('depth', required_type='float')
            still_water_pos = field_inp.get_value('still_water_position', required_type='float')

            # Project the colour function to DG0 (set degree to -1 to prevent this)
            # The special projection in BaseWaveField handles partial filling of cells
            # by use of quadrature elements in the C++ expression and a DG0 projection
            self.colour_projection_degree = 6
            self.colour_projection_form = None

            # Define the C++ code (YOU MUST CHANGE THESE TO SOME VALID C++ EXPRESSIONS)
            # Use x[0] for the horizontal coordinate and x[2] for the vertical coordinate,
            # conversion to 2D is done below if the simulation is 2D
            cpp_e = 'C++ code for the wave elevation;'
            cpp_u = 'C++ code for the particle velocity in the horizontal direction;'
            cpp_w = 'C++ code for the particle velocity in the vertical direction;'

            # Store the C++ code
            self._cpp['elevation'] = cpp_e
            self._cpp['uhoriz'] = cpp_u
            self._cpp['uvert'] = cpp_w
            self._cpp['c'] = 'x[2] <= (%s) ? 1.0 : 0.0' % cpp_e

            # Adjust the z-coordinate such that the bottom is at z=0
            # (make changes if your code assumes that the free surface is at z=0)
            zdiff = still_water_pos - self.h
            for k, v in list(self._cpp.items()):
                self._cpp[k] = v.replace('x[2]', '(x[2] - %r)' % zdiff)

            # Adjust the C++ code z-coordinate if the simulation is 2D (then z -> y)
            if self.simulation.ndim == 2:
                for k, v in list(self._cpp.items()):
                    self._cpp[k] = v.replace('x[2]', 'x[1]')


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
