Programmers guide
=================

This documentation shows the API of the main classes and gives information 
about scripting Ocellaris. You can also find help to understand the main
code base if you want to help further develop Ocellaris.   

.. toctree::

    simulation
    boundary_conditions
    solver
    scripting_and_interactive_console

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
line by running the :meth:`ocellaris.simulation.Input.input.read_yaml` method.
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

A simplified replification of the above in a script would be:

.. code-block:: python

    from ocellaris import Simulation, run_simulation

    sim = Simulation()
    sim.input.read_yaml('template.inp')
    run_simulation(sim)

Read more about scripting in the scripting section
