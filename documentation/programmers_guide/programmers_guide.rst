Programmers guide
=================

This page contains some pointers on how to understand the Ocellaris code.   

.. contents:: :local:


A brief introduction
--------------------

The following is a description of what happens when the user starts Ocellaris
by running the following on the command line::

    python -m ocellaris INPUT_FILE
    
Ocellaris starts by running the :func:`main()` function in the 
:mod:`ocellaris.__main__` module. This function will create an object of the
:class:`ocellaris.Simulation` class. This simulation object will be central to
the execution of Ocellaris and it will be passed around to allmost all pieces
of the code. Everyone who wants to look at the input or access the calculated
solution must do this through the simulation class.

The main function will now read the input file given by the user on the command
line by running the :meth:`ocellaris.simulation.Input.input.read_json` method.
The code will also set up logging / console output and print a banner unless
the user has set the log level so high that INFO messages will not be printed.

Next the :func:`ocellaris.run_simulation` function is called and then the 
:mod:`ocellaris.__main__` module will take no more part in the running of 
Ocellaris except for printing a goodbye message at the end.

The main task of setting up and running the simulation is done in the
:mod:`ocellaris.run` module. This is where the :func:`ocellaris.run_simulation`
function is implemented along with several utility functions. The following
actions are performed here:

- Load the mesh
- Create function spaces
- Create boundary conditions
- Load physical constants
- Create the multiphase model (controls density and viscosity)
- Create probes which can report solution data to file and/or show interactive 
    plots during the simulation
- Populate the :attr:`ocellaris.Simulation.data` dictionary with the mesh,
  function spaces, boundary contitions etc
- Create the solver
- Run the solver
- Report how long each part of the simulation took


Simulation classes
------------------

The simulation classes holds the simulation data (velocity, pressure, function 
spaces, mesh etc) and is also responsible for most utility functionality such
as plugins (hooks), logging, reporting, plotting and input file handling. 

.. autoclass:: ocellaris.Simulation

    .. attribute:: input
        :annotation:
        
        An :class:`ocellaris.simulation.Input` object that holds the input from
        the input file provided by the user
    
    .. attribute:: hooks
        :annotation:
        
        An :class:`ocellaris.simulation.Hooks` object that keeps track of
        functions that should run at certain times during the simulation

    .. attribute:: plotting
        :annotation:
        
        An :class:`ocellaris.simulation.Plotting` object that helps creating 
        plots of the solution while the simulation is running

    .. attribute:: reporting
        :annotation:
        
        An :class:`ocellaris.simulation.Reporting` object that helps report
        summary values each time step

    .. attribute:: log
        :annotation:
        
        An :class:`ocellaris.simulation.Log` object that helps with logging
        messages to screen while the simulation is running

.. autoclass:: ocellaris.simulation.Input

.. autoclass:: ocellaris.simulation.Hooks

.. autoclass:: ocellaris.simulation.Plotting

.. autoclass:: ocellaris.simulation.Reporting

.. autoclass:: ocellaris.simulation.Log


Boundary conditions
-------------------

The boundary condition code will both identify regions of the boundary given by
the user in the input file and create boundary condition objects for each 
function (velocity, pressure ...) in this region.

.. autoclass:: ocellaris.boundary_conditions.BoundaryRegion

.. autoclass:: ocellaris.boundary_conditions.BoundaryCondition


The solver
----------

The solver uses the simulation classes and runs a time loop to solve the time
dependent Navier-Stokes equations.

.. autoclass:: ocellaris.solvers.ipcs.SolverIPCS


.. _sec-interactive-console:

Interactive console
-------------------

At the end of each time step Ocellaris will optionally open an interactive
console so that you can inspect the internal state of the simulation. To
access this pres :kbd:`d` then :kbd:`Enter` ("d" for debug). At the end of the
next time step the console should open and you will have full access to the
internal variables. The variables are listed so that you can get a head start

If you press :kbd:`Ctrl+d` inside the interactive console Ocellaris will
continue running the time loop. If you type ``exit()`` or  ``quit()`` you will
stop Ocellaris and return to the command line immediately. 

It is also possible to specify that the console should open at the end of the
simulation. If you want this put the following on the input file:

.. code-block:: javascript

    "console_at_end": true

This can be very useful for ad-hoc postprocessing of the simulation results.