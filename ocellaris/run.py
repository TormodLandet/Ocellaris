import sys, traceback, time
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
        
        simulation.setup()
        success = True
    
    except OcellarisError as e:
        simulation.log.error('ERROR === '*8)
        simulation.log.error('\n%s\n\n%s\n' % (e.header, e.description))
    except KeyboardInterrupt as e:
        simulation.log.error('========== You pressed Ctrl+C -- STOPPING ==========')
    except BaseException as e:
        simulation.log.error('=== EXCEPTION =='*5)    
        tb = traceback.format_tb(sys.exc_info()[2])
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb))
        e_type = type(e).__name__
        simulation.log.error('Got %s exception when running setup:\n%s' % (e_type, str(e)))
    
    # Check if the setup ran without problems
    if not success and not catch_exceptions:
        raise # Re-raise the exception gotten from running the solver
    
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
    simulation.log.info("\nRunning simulation on %d CPUs...\n" % simulation.ncpu)
    simulation.t_start = time.time()
    
    try:
        success = False
        simulation.solver.run()
        success = True
    except OcellarisError as e:
        simulation.hooks.simulation_ended(success)
        simulation.log.error('ERROR === '*8)
        simulation.log.error('\n%s\n\n%s\n' % (e.header, e.description))
    except KeyboardInterrupt as e:
        simulation.hooks.simulation_ended(success)
        simulation.log.error('========== You pressed Ctrl+C -- STOPPING ==========')
    except SystemExit as e:
        simulation.success = False  # this is just used for debugging, no fancy summary needed
        simulation.log.error('========== SystemExit - exit() was called ==========')
    except BaseException as e:
        simulation.hooks.simulation_ended(success)
        simulation.log.error('=== EXCEPTION =='*5)    
        tb = traceback.format_tb(sys.exc_info()[2])
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb))
        e_type = type(e).__name__
        simulation.log.error('Got %s exception when running solver:\n%s' % (e_type, str(e)))
    
    # Check if the solver ran without problems
    if not success and not catch_exceptions:
        raise # Re-raise the exception gotten from running the solver 
    
    ##############################################################################################
    # Limited support for postprocessing implemented below. It is generally better to use Paraview
    # or similar tools on the result files instead of using the dolfin plot commands here 
        
    # Show dolfin plots?
    if simulation.input.get_value('output/plot_at_end', False, 'bool'):
        plot_at_end(simulation)
    
    # Optionally show the console for debugging and ad-hoc posprocessing
    console_at_end = simulation.input.get_value('console_at_end', False, 'bool')
    console_on_error = simulation.input.get_value('console_on_error', False, 'bool')
    if console_at_end  or (not success and console_on_error):
        run_debug_console(simulation)
        
    return success


def plot_at_end(simulation):
    """
    Open dolfin plotting windows with results at the end of
    the simulation
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
    
    dolfin.interactive()
