import sys
import code

def debug_console_hook(simulation):
    """
    Start the debug console if the user writes "d"+Enter
    """
    # The select() system call does not work on windows
    if 'win' in sys.platform:
        return
    
    # Check if there is input on stdin. If there is a line
    # containing only a "d" then start the debuc console
    import select
    has_input = lambda: sys.stdin in select.select([sys.stdin], [], [], 0)[0]
    while has_input(): 
        line = sys.stdin.readline()
        if line.strip().lower() == 'd':
            return run_debug_console(simulation)

def run_debug_console(simulation):
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
    for name in sorted(names):
        banner.append(' - %s' % name)
    banner.append('\nPress Ctrl+D to continue running the simulation.'
                  '\nRunning exit() or quit() will stop Ocellaris.')
    banner.append('\n\n>>> from dolfin import *')
    
    print '=== OCELLARIS CONSOLE === '*3
    code.interact('\n'.join(banner), local=debug_locals)
    print '= OCELLARIS CONSOLE ===='*3
