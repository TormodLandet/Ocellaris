.. _scripting-ocellaris:

Scripting Ocellaris
===================

There are three main ways of scripting Ocellaris:

#) Generate input files by a script and run Ocellaris on these
#) Include Python code in the input file (see the :ref:`User Guide
   <user_guide>` for details)
#) Use the Python programming interface

A brief example of the first would be something like this that runs Ocellaris
with two different time steps by producing input files and starting Ocellaris
as a command line application:

.. code-block:: python

    import yaml, subprocess

    # Load the template input file
    with open('template.inp', 'rt') as inp:
        input_dict = yaml.load(inp)

    for dt in [0.1, 0.05]:
        # Modify the input file
        prefix = 'test_dt_%.3f' % dt
        input_dict['time']['dt'] = dt
        input_dict['output']['prefix'] = prefix

        # Save the modified input file
        new_inp_file = prefix + '.inp'
        with open(new_inp_file, 'wt') as out:
            yaml.dump(input_dict, out)

        # Run Ocellaris with the modified input file
        p = subprocess.Popen(['python', '-m', 'ocellaris', new_inp_file],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()


The same can be accomplished from a script which uses the Ocellaris Python API:

.. code-block:: python

    from ocellaris import Simulation, setup_simulation, run_simulation

    for dt in [0.1, 0.05]:
        prefix = 'test_dt_%.3f' % dt

        # Create a simulation object, load an input file and modify the time step
        sim = Simulation()
        sim.input.read_yaml('template.inp')
        sim.input.set_value('time/dt', dt)
        sim.input.set_value('output/prefix', prefix)

        # Run Ocellaris
        setup_simulation(sim)
        run_simulation(sim)

For more information about what you can do with the simulation object, see the
:ref:`simulation-api` documentation.

Examples of this can be seen in the convergence scripts which can be found in
the ``cases/`` subdirectory of the Ocellaris repository.
