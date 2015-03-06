import dolfin
from ocellaris import get_version, get_detailed_version, Simulation, run_simulation

# Names for the available log levels in Dolfin and Ocellaris
AVAILABLE_LOG_LEVELS = {'critical': dolfin.CRITICAL,
                        'error': dolfin.ERROR,
                        'warning': dolfin.WARNING,
                        'info': dolfin.INFO,
                        'progress': dolfin.PROGRESS,
                        'debug': dolfin.DEBUG}

def main(inputfile):
    """
    Run Ocellaris
    """
    sim = Simulation()
    sim.input.read_json(inputfile)
    sim.log.setup()
    
    # Set the Ocellaris log level
    log_level = sim.input.get_value('output/ocellaris_log_level', 'info')
    sim.log.set_log_level(AVAILABLE_LOG_LEVELS[log_level])
    
    # Set the Dolfin log level
    df_log_level = sim.input.get_value('output/dolfin_log_level', 'warning')
    dolfin.set_log_level(AVAILABLE_LOG_LEVELS[df_log_level])
    
    version = get_detailed_version() or get_version()
    
    sim.log.info('='*80)
    sim.log.info('                  Ocellaris   %s' % version) 
    sim.log.info('='*80)
    sim.log.info()
    
    run_simulation(sim)
    
    sim.log.info('='*80)
    sim.log.info('Ocellaris finished successfully')

if __name__ == '__main__':
    # Get command line arguments
    import argparse
    parser = argparse.ArgumentParser(prog='ocellaris',
                                     description='Discontinuous Galerkin Navier-Stokes solver')
    parser.add_argument('inputfile', help='Name of file containing simulation '
                        'configuration on the Ocellaris json input format')
    args = parser.parse_args()
    
    # Run Ocellaris
    main(args.inputfile)
