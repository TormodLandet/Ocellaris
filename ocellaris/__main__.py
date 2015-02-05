import time
from ocellaris import version, Simulation, run_simulation

def main(inputfile):
    """
    Run Ocellaris
    """
    t1 = time.time()
    print 'Ocellaris v %s' % version 
    print '========================================'
    print
    
    sim = Simulation()
    sim.read_json_input_file(inputfile)
    run_simulation(sim)
    
    print '========================================'
    print 'Ocellaris DONE after %.3f seconds' % (time.time() - t1)

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
