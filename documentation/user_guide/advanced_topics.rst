Advanced topics
===============


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
        description = "This is mu custom solver!"

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
