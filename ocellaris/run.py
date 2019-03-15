# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import sys
import os
import traceback
import time
from .utils import OcellarisError, run_debug_console
import dolfin


def setup_simulation(simulation, setup_logging=True, catch_exceptions=False):
    """
    Prepare and run a simulation

    This will run the "real" code ``run_simulation_without_error_handling``
    in a way that handles errors in a user friendly manner (log them etc)

    :type simulation: ocellaris.Simulation
    """
    try:
        success = False

        if setup_logging:
            simulation.log.setup()

        # Validate the input file format - try to catch any misspellings
        validate_input_file(simulation)

        # Run the setup routines - can take a while on large meshes
        simulation.setup()
        success = True

    except OcellarisError as e:
        simulation.log.error('ERROR === ' * 8)
        simulation.log.error('\n%s\n\n%s\n' % (e.header, e.description))
        tb, err = sys.exc_info()[2], e
    except KeyboardInterrupt as e:
        simulation.log.error('========== You pressed Ctrl+C -- STOPPING ==========')
        tb, err = sys.exc_info()[2], e
        tb_msg = traceback.format_tb(tb)
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb_msg))
        simulation.log.error(
            'You pressed Ctrl+C or setup got a SIGINT/SIGTERM∕SIGQUIT signal'
        )
    except BaseException as e:
        simulation.log.error('=== EXCEPTION ==' * 5)
        tb, err = sys.exc_info()[2], e
        tb_msg = traceback.format_tb(tb)
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb_msg))
        e_type = type(e).__name__
        simulation.log.error(
            'Got %s exception when running setup:\n%s' % (e_type, str(e))
        )

    # Check if the setup ran without problems
    if not success and not catch_exceptions:
        raise err.with_traceback(
            tb
        )  # Re-raise the exception gotten from running the setup

    return success


def run_simulation(simulation, catch_exceptions=False):
    """
    Prepare and run a simulation

    This will run the "real" code ``run_simulation_without_error_handling``
    in a way that handles errors in a user friendly manner (log them etc)

    :type simulation: ocellaris.Simulation
    """
    # Print information about configuration parameters
    simulation.log.info('\nSimulation configuration when starting the solver:')
    simulation.log.info('\n{:-^80}'.format(' configuration begin '))

    simulation.log.info(str(simulation.input))
    simulation.log.info('{:-^80}'.format(' configuration end '))
    simulation.log.info("\nCurrent time: %s" % time.strftime('%Y-%m-%d %H:%M:%S'))
    simulation.log.info(
        "Degrees of freedom: %d" % simulation.ndofs
        + (
            ''
            if simulation.ncpu == 1
            else ' (%d per process)' % (simulation.ndofs / simulation.ncpu)
        )
    )
    simulation.log.info("\nRunning simulation on %d CPUs...\n" % simulation.ncpu)
    simulation.t_start = time.time()

    try:
        success = False
        simulation.flush()
        simulation.solver.run()
        success = True
    except OcellarisError as e:
        simulation.log.error('ERROR === ' * 8)
        simulation.log.error('\n%s\n\n%s\n' % (e.header, e.description))
        simulation.hooks.simulation_ended(success)
        tb, err = sys.exc_info()[2], e
    except KeyboardInterrupt as e:
        # This will also handle signals such as SIGTERM and SIGQUIT in addition to SIGINT
        # due to the signal handler setup in __main___.py
        simulation.log.error('========== You pressed Ctrl+C -- STOPPING ==========')
        tb, err = sys.exc_info()[2], e
        tb_msg = traceback.format_tb(tb)
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb_msg))
        simulation.log.error(
            'You pressed Ctrl+C or the solver got a SIGINT/SIGTERM∕SIGQUIT signal'
        )
        if simulation.input.get_value('output/save_restart_file_at_end', True, 'bool'):
            simulation.io.last_savepoint_is_checkpoint = True
            simulation.hooks.simulation_ended(success)
    except SystemExit as e:
        simulation.success = (
            False
        )  # this is just used for debugging, no fancy summary needed
        simulation.log.error('========== SystemExit - exit() was called ==========')
        tb, err = sys.exc_info()[2], e
    except BaseException as e:
        simulation.log.error('=== EXCEPTION ==' * 5)
        tb, err = sys.exc_info()[2], e
        tb_msg = traceback.format_tb(tb)
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb_msg))
        e_type = type(e).__name__
        simulation.log.error(
            'Got %s exception when running solver:\n%s' % (e_type, str(e))
        )
        simulation.hooks.simulation_ended(success)
    simulation.flush()

    # Check if the solver ran without problems
    if not success and not catch_exceptions:
        raise err.with_traceback(
            tb
        )  # Re-raise the exception gotten from running the solver

    ##############################################################################################
    # Limited support for postprocessing implemented below. It is generally better to use Paraview
    # or similar tools on the result files instead of using the dolfin plot commands here

    # Show dolfin plots?
    if simulation.input.get_value('output/plot_at_end', False, 'bool'):
        plot_at_end(simulation)

    # Optionally show the console for debugging and ad-hoc posprocessing
    console_at_end = simulation.input.get_value('console_at_end', False, 'bool')
    console_on_error = simulation.input.get_value('console_on_error', False, 'bool')
    if console_at_end or (not success and console_on_error):
        run_debug_console(simulation)

    return success


def validate_input_file(simulation):
    """
    Use YSchema to validate the input file
    """
    simulation.log.info('Validating input using YSchema')
    try:
        import yschema
    except ImportError:
        simulation.log.warning('Could not import the "yschema" Python library')
        return

    mydir = os.path.abspath(os.path.dirname(__file__))
    schemafile = os.path.join(mydir, 'input_file_schema.yml')
    with open(schemafile, 'rt') as yml:
        schema = yschema.yaml_ordered_load(yml)

    # Run the validation routines
    errors = yschema.validate(simulation.input, schema, return_errors=True)
    if not errors:
        simulation.log.info('    Validation status OK')
        return

    simulation.log.error(
        '    The input file is not valid! Ocellaris will still'
        ' continue since this may be a false negative (the'
        ' validation schema is not perfect), but beware!'
    )
    prefix1 = '    ERROR: '
    prefix2 = ' ' * len(prefix1)
    for error in errors:
        simulation.log.warning(prefix1 + error.replace('\n', '\n' + prefix2))


def plot_at_end(simulation):
    """
    Open dolfin plotting windows with results at the end of
    the simulation. Only useful for simple 2D cases
    """
    # Plot velocity components
    for d in range(simulation.ndim):
        name = 'u%d' % d
        dolfin.plot(simulation.data[name], title=name)
        name = 'u_star%d' % d
        dolfin.plot(simulation.data[name], title=name)

    dolfin.plot(simulation.data['u'], title='u')

    # Plot pressure
    dolfin.plot(simulation.data['p'], title='p')

    # Plot colour function
    if 'c' in simulation.data:
        dolfin.plot(simulation.data['c'], title='c')

    dolfin.plot(simulation.data['boundary_marker'], title='boundary_marker')

    from matplotlib import pyplot

    pyplot.show()
