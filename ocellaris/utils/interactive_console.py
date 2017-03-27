import sys
import code
import readline
import rlcompleter
import cProfile, pstats
from .timer import log_timings


def interactive_console_hook(simulation):
    """
    Start the interactive console if the user writes a valid
    command to the prompt while Ocellaris is running its time
    loop, ie. writing "[d] + [Enter]" will start the debug
    console.
    
    Commands:
    
    * d - debug console
    * p - plot fields
    * r - write restart file
    * s - stop the simulation (changes the maximum simulation
          time to current time)
    * t - print timings to screen
    * prof N - profile the next N time steps
    
    Interactive console commands are not available on Windows
    or during non-interactive (queue system/batch) use
    """
    # Check if the user has written something on stdin for us
    if simulation.rank == 0:
        commands = get_input_from_terminal()
    else:
        commands = []
    
    # Make sure all processes get the same commands
    if simulation.ncpu > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        commands = comm.bcast(commands)
    
    for command in commands:
        if command == 'd' and simulation.ncpu == 1:
            # d == "debug" -> start the debug console
            return run_debug_console(simulation)
        
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
            simulation.log.info('\nCommand line action:\n  Setting simulation '
                                'control parameter tmax to %r\n' % simulation.time)
            simulation.input['time']['tmax'] = simulation.time

        elif command == 't':
            # t == "timings" -> show timings
            simulation.log.info('\nCommand line action:\n  Showing timings')
            log_timings(simulation)
        
        elif command.startswith('prof') and simulation.rank == 0:
            # Run the profiler
            try:
                num_timesteps = int(command.split()[1])
            except:
                simulation.log.warning('Did not understand requested number of profile time steps')
                return
            simulation._profile_after_n_timesteps = num_timesteps+1
            simulation._profile_object = cProfile.Profile()
            simulation._profile_object.enable()
            simulation.log.info('\nCommand line action:\n  Starting profile')
    
    # Write the profile information to file after N time steps
    if hasattr(simulation, '_profile_after_n_timesteps'):
        simulation._profile_after_n_timesteps -= 1
        simulation.log.info('Profile will end after %d time steps' % simulation._profile_after_n_timesteps)
        if simulation._profile_after_n_timesteps == 0:
            simulation._profile_object.disable()
            stats = pstats.Stats(simulation._profile_object)
            stats.strip_dirs()
            stats.sort_stats('cumulative')
            stats.print_stats(30)
            print 'Saving cProfile trace to "prof.out"'
            stats.dump_stats('prof.out')
            del simulation._profile_after_n_timesteps
            del simulation._profile_object


def get_input_from_terminal():
    """
    Read stdin to see if there are some commands for us to execute
    """
    # The select() system call does not work on windows
    if 'win' in sys.platform:
        return []
    
    # If we are not running interactively there will be no commands
    if not sys.__stdin__.isatty():
        return []
    
    import select
    
    commands = []
    has_input = lambda: sys.stdin in select.select([sys.stdin], [], [], 0)[0]
    
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
    names = simulation.data.keys() + ['simulation']
    for i, name in enumerate(sorted(names)):
        if i % 4 == 0:
            banner[-1] += '\n'
        banner[-1] += '  %-18s' % name
    banner.append('\n\nPress Ctrl+D to continue running the simulation.'
                  '\nRunning exit() or quit() will stop Ocellaris.')
    
    # Add some convenience functions
    funcs, info = define_convenience_functions(simulation)
    debug_locals.update(funcs)
    banner.extend(info)
    
    # Setup tab completion
    readline.set_completer(rlcompleter.Completer(debug_locals).complete)
    readline.parse_and_bind("tab: complete")
    
    if not show_banner:
        banner = []
    
    print '=== OCELLARIS CONSOLE === '*3
    banner.append('\n>>> from dolfin import *')
    code.interact('\n'.join(banner), local=debug_locals)
    print '= OCELLARIS CONSOLE ===='*3


def define_convenience_functions(simulation):
    """
    Some functions that are nice to have for debugging purposes
    """
    import dolfin
    
    info = []
    funcs = {}
    
    # Convenience plotting function
    fields = [name for name in ('u', 'p', 'c', 'p_hydrostatic') if name in simulation.data]
    info.append('Running plot_all() will plot %s' % ' & '.join(fields))
    def plot_all():
        for name in fields:
            field = simulation.data[name]
            dolfin.plot(field, title=name, tag=name)
        dolfin.interactive()
    funcs['plot_all'] = plot_all
    
    return funcs, info
