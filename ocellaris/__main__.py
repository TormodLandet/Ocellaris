# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import sys
import os
from ocellaris import get_detailed_version, Simulation, setup_simulation, run_simulation
import platform
import dolfin
import meshio
import yaml
import mpi4py
import h5py
from petsc4py import PETSc, __version__ as petsc4py_version


# Restore signals in non-interactive background shells
import signal

signal.signal(signal.SIGINT, signal.default_int_handler)
signal.signal(signal.SIGTERM, signal.default_int_handler)
signal.signal(signal.SIGQUIT, signal.default_int_handler)


def main(inputfile, input_override):
    """
    Run Ocellaris
    """
    if os.environ.get('OCELLARIS_SUPER_DEBUG', False):
        print('FOUND OCELLARIS_SUPER_DEBUG in environment')
        from ocellaris.utils.debug import enable_super_debug

        enable_super_debug()

    sim = Simulation()

    # Read input
    if sim.io.is_restart_file(inputfile):
        sim.io.load_restart_file_input(inputfile)
    else:
        sim.input.read_yaml(inputfile)

    # Alter input by values given on the command line
    override_input_variables(sim, input_override)

    # Setup logging before we start printing anything
    sim.log.setup()

    # Print banner with Ocellaris version number
    version = get_detailed_version()
    location = os.path.split(os.path.abspath(__file__))[0]
    sim.log.info('=' * 80)
    sim.log.info('                  Ocellaris   %s' % version)
    sim.log.info('=' * 80)
    sim.log.info('Installed at:')
    sim.log.info('    %s' % location)
    sim.log.info('    host: %s' % platform.node())
    sim.log.info()

    # Print some version information
    sim.log.info('Running on Python %s' % sys.version)
    sim.log.info('    Using dolfin %s' % dolfin.__version__)
    sim.log.info('    Using mpi4py %s' % mpi4py.__version__)
    sim.log.info('    Using h5py   %s' % h5py.__version__)
    sim.log.info('    Using meshio %s' % meshio.__version__)
    sim.log.info('    Using PyYAML %s' % yaml.__version__)
    sim.log.info('    Using petsc4py %s' % petsc4py_version)
    sim.log.info('    Using PETSc %d.%d.%d\n' % PETSc.Sys.getVersion())

    # Setup the Ocellaris simulation
    ok = setup_simulation(sim, setup_logging=False, catch_exceptions=True)
    if not ok:
        sim.log.error('Setup did not suceed, exiting')
        sys.exit(1)

    if sim.restarted:
        # Load previous results
        sim.io.load_restart_file_results(inputfile)

    # Run the Ocellaris simulation time loop
    run_simulation(sim, catch_exceptions=True)

    sim.log.info('=' * 80)
    if sim.success:
        sim.log.info('Ocellaris finished successfully')
    else:
        sim.log.info('Ocellaris finished with errors')


def override_input_variables(simulation, input_override):
    """
    The user can override values given on the input file via
    command line parameters like::

        --set-input time/dt=0.1

    This code updates the input dictionary with these modifications
    """
    if input_override is None:
        return

    for overrider in input_override:
        # Overrider is something like "time/dt=0.1"
        path = overrider.split('=')[0]
        value = overrider[len(path) + 1 :]

        # Find the correct input sub-dictionary
        path_elements = path.split('/')
        base_path = '/'.join(path_elements[:-1])
        base = simulation.input.ensure_path(base_path)

        # Convert value to Python object
        try:
            py_value = eval(value)
        except Exception as e:
            print('ERROR: Input variable given via command line argument failed:')
            print('ERROR:       --set-input "%s"' % overrider)
            print('ERROR: Got exception: %s' % str(e))
            exit(-1)

        # The last path element may be a list index
        try:
            idx = int(path_elements[-1])
        except ValueError:
            idx = path_elements[-1]

        # Update the input sub-dictionary
        base[idx] = py_value
        simulation.log.info('Command line overriding %r in %s to %r' % (idx, base_path, py_value))


def run_from_console():
    # Get command line arguments
    import argparse

    parser = argparse.ArgumentParser(
        prog='ocellaris', description='Discontinuous Galerkin Navier-Stokes solver'
    )
    parser.add_argument(
        'inputfile',
        help='Name of file containing simulation '
        'configuration on the Ocellaris YAML input format',
    )
    parser.add_argument(
        '--set-input',
        '-i',
        action='append',
        help='Set an input key. Can be added several '
        'times to set multiple input keys. Example: --set-input time/dt=0.1',
    )

    parser.add_argument(
        '--pystuck',
        action='store_true',
        help='Activate pystuck to debug hangs (only on MPI rank 0)',
    )

    args = parser.parse_args()

    # Enable debuging of stuck processes
    if args.pystuck and dolfin.MPI.comm_world.rank == 0:
        try:
            import pystuck

            pystuck.run_server()
        except Exception as e:
            print('Could not start pystuck')
            print(e)
            print('Starting Ocellaris without pystuck')

    # Run Ocellaris
    main(args.inputfile, args.set_input)


if __name__ == '__main__':
    run_from_console()
