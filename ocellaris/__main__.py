import time
from ocellaris import get_version, get_detailed_version, Simulation, run_simulation

def main(inputfile):
    """
    Run Ocellaris
    """
    version = get_detailed_version() or get_version()
    
    print '='*80
    print '                  Ocellaris   %s' % version 
    print '='*80
    print
    
    sim = Simulation()
    sim.input.read_json(inputfile)
    run_simulation(sim)
    
    print '='*80
    print 'Ocellaris exiting successfully'

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
