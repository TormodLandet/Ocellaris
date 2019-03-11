# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os
import sys
import code
import readline
import rlcompleter
import cProfile
import pstats
from .timer import log_timings


def interactive_console_hook(simulation):
    """
    Start the interactive console if the user writes a valid
    command to the prompt while Ocellaris is running its time
    loop, ie. writing "[d] + [Enter]" will start the debug
    console.

    Commands:

    * d - debug console
    * f - flush open files
    * i - change input
    * p - plot fields
    * prof N - profile the next N time steps
    * r - write restart file
    * s - stop the simulation (changes the maximum simulation
          time to current time)
    * t - print timings to screen
    * w - write output file (vtk or xdmf)

    Interactive console commands are not available on Windows
    or during non-interactive (queue system/batch) use
    """
    # Command file - used in addition to stdin if it exits
    command_file = simulation.input.get_output_file_path(
        'command_control_file', '.COMMANDS'
    )

    # Check if the user has written something on stdin for us
    if simulation.rank == 0:
        commands = read_stdin_and_cmdfile(command_file)
    else:
        commands = []

    # Make sure all processes get the same commands
    if simulation.ncpu > 1:
        comm = simulation.data['mesh'].mpi_comm()
        commands = comm.bcast(commands)

    for command in commands:
        cwords = command.split()

        # One-word commands in alphabetical order:

        if command == 'd' and simulation.ncpu == 1:
            # d == "debug" -> start the debug console
            run_debug_console(simulation)

        elif command == 'f':
            # f == "flush" -> flush open files
            simulation.log.info('\nCommand line action:\n  Flushing open files')
            simulation.hooks.run_custom_hook('flush')

        elif command == 'p':
            # p == "plot" -> plot field variable
            funcs, _ = define_convenience_functions(simulation)
            simulation.log.info('\nCommand line action:\n  Plotting fields')
            funcs['plot_all']()

        elif command == 'r':
            # r == "restart" -> write restart file
            simulation.log.info('\nCommand line action:\n  Writing restart file')
            simulation.io.write_restart_file()

        elif command == 's':
            # s == "stop" -> stop the simulation
            simulation.log.info(
                '\nCommand line action:\n  Setting simulation '
                'control parameter tmax to %r\n' % simulation.time
            )
            simulation.input['time']['tmax'] = simulation.time

        elif command == 't':
            # t == "timings" -> show timings
            simulation.log.info('\nCommand line action:\n  Showing timings')
            log_timings(simulation)

        # Multi-word commands in alphabetical order:

        elif len(cwords) < 2:
            continue

        elif cwords[0] == 'i':
            # i == "input" -> set input variable, e.g., "i time/dt = 0.04"
            assignment = ' '.join(cwords[1:])
            simulation.log.info(
                '\nCommand line action:\n  Input modification: %s' % assignment
            )
            i = assignment.find('=')
            if i == -1:
                simulation.log.warning('Malformed input, no "=" found')
                continue
            path = assignment[:i].strip()
            value = assignment[i + 1 :].strip()
            try:
                py_value = eval(value)
            except Exception:
                simulation.log.warning('Could not get Python value from %r' % value)
                continue
            simulation.log.info('  Setting input %r to %r' % (path, py_value))
            simulation.input.set_value(path, py_value)

        elif cwords[0] == 'w':
            # w == "write" -> write output file (restart files use "r" command)
            file_format = cwords[1]
            simulation.log.info(
                '\nCommand line action:\n  Writing %s file' % file_format
            )
            if file_format == 'vtk':
                simulation.io.lvtk.write()
            elif file_format == 'xdmf':
                simulation.io.xdmf.write()
            else:
                simulation.log.warning('Unsupported file type %r' % file_format)

        elif cwords[0] == 'prof' and simulation.rank == 0:
            # Run the profiler
            try:
                num_timesteps = int(cwords[1])
            except Exception:
                simulation.log.warning(
                    'Did not understand requested number of profile time steps'
                )
                return
            simulation._profile_after_n_timesteps = num_timesteps + 1
            simulation._profile_object = cProfile.Profile()
            simulation._profile_object.enable()
            simulation.log.info('\nCommand line action:\n  Starting profile')

    # Write the profile information to file after N time steps
    if hasattr(simulation, '_profile_after_n_timesteps'):
        simulation._profile_after_n_timesteps -= 1
        simulation.log.info(
            'Profile will end after %d time steps'
            % simulation._profile_after_n_timesteps
        )
        if simulation._profile_after_n_timesteps == 0:
            simulation._profile_object.disable()
            stats = pstats.Stats(simulation._profile_object)
            stats.strip_dirs()
            stats.sort_stats('cumulative')
            stats.print_stats(30)
            print('Saving cProfile trace to "prof.out"')
            stats.dump_stats('prof.out')
            del simulation._profile_after_n_timesteps
            del simulation._profile_object


def read_stdin_and_cmdfile(command_file):
    """
    Read stdin to see if there are some commands for us to execute.
    We also check the file output_prefix.COMMANDS and read this
    if it exists for cases where stdin is not easily accessible
    for the user. The command file is DELETED after reading
    """
    # Read commands from the control command file if it exists
    # and then delete it after reading to avoid re-running the
    # commands every time step
    commands = []
    if os.path.isfile(command_file):
        try:
            with open(command_file) as inp:
                lines = inp.readlines()
            os.remove(command_file)
            commands.extend(cmd.strip() for cmd in lines)
        except Exception:
            pass

    # We need the select() system call to check if is possible
    # to read from stdin without blocking (waiting for input).
    # The select() system call does not work on windows
    if 'win' in sys.platform:
        return commands

    # If we are not running interactively there will be no commands
    if not sys.__stdin__.isatty():
        return commands

    # Check if we can read from stdin without blocking
    import select

    def has_input():
        return sys.stdin in select.select([sys.stdin], [], [], 0)[0]

    # Check if there is input on stdin and read it if it exists
    while has_input():
        line = sys.stdin.readline()
        command = line.strip().lower()
        commands.append(command)

    return commands


def run_debug_console(simulation, show_banner=True):
    """
    Run a debug console with some useful variables available
    """
    banner = ['Ocellaris interactive console\n']

    # Everything from dolfin should be available
    import dolfin

    debug_locals = dict(**dolfin.__dict__)

    # All variables in the data dict should be available
    debug_locals.update(simulation.data)

    # The simulation object shoule be available
    debug_locals['simulation'] = simulation

    # Create a banner to show before the console
    banner.append('Available variables:')
    names = list(simulation.data.keys()) + ['simulation']
    for i, name in enumerate(sorted(names)):
        if i % 4 == 0:
            banner[-1] += '\n'
        banner[-1] += '  %-18s' % name
    banner.append(
        '\n\nPress Ctrl+D to continue running the simulation.'
        '\nRunning exit() or quit() will stop Ocellaris.'
    )

    # Add some convenience functions
    funcs, info = define_convenience_functions(simulation)
    debug_locals.update(funcs)
    banner.extend(info)

    # Setup tab completion
    readline.set_completer(rlcompleter.Completer(debug_locals).complete)
    readline.parse_and_bind("tab: complete")

    if not show_banner:
        banner = []

    print('=== OCELLARIS CONSOLE === ' * 3)
    banner.append('\n>>> from dolfin import *')
    code.interact('\n'.join(banner), local=debug_locals)
    print('= OCELLARIS CONSOLE ====' * 3)


def define_convenience_functions(simulation):
    """
    Some functions that are nice to have for debugging purposes
    """
    import dolfin
    from matplotlib import pyplot

    info = []
    funcs = {}

    # Convenience plotting function
    fields = [
        name for name in ('u', 'p', 'c', 'p_hydrostatic') if name in simulation.data
    ]
    info.append('Running plot_all() will plot %s' % ' & '.join(fields))

    def plot_all():
        for name in fields:
            field = simulation.data[name]
            dolfin.plot(field, title=name, tag=name)
        pyplot.show()

    funcs['plot_all'] = plot_all
    funcs['pyplot'] = pyplot

    return funcs, info
