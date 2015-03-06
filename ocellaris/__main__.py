from ocellaris import get_version, get_detailed_version, Simulation, run_simulation

def main(inputfile):
    """
    Run Ocellaris
    """
    sim = Simulation()
    sim.input.read_json(inputfile)
    sim.log.setup()
    
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
