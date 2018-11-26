Utility scripts
===============

Ocellaris comes with utility scripts for automating common tasks in the
``scripts`` directory.


.. contents:: Contents
    :local:


orun.py - Run Ocellaris on a HPC with automatic restarts
--------------------------------------------------------

For running Ocellaris on a cluster. MPI solvers can be fragile and sometimes
will hang for no apparent reason. If you run your simulations with orun.py and
make sure to write restart files every once in a while then orun will monitor
your output and kill and restart the simulation if it appears to have stopped.

An example SLURM job script using orun.py:

.. code-block:: bash

    #!/bin/bash

    # Metadata:
    #SBATCH --job-name=Ocellaris+FEniCS
    #SBATCH --account=XXXXXXXXX
    #SBATCH --output="slurm-%j.STDOUT"
    #SBATCH --error="slurm-%j.STDERR"
    
    # Resource requests:
    #SBATCH --time=3-00:00:00
    #SBATCH --mem-per-cpu=3936
    #SBATCH --ntasks-per-node=16 --nodes=2

    # Abel cluster specific
    set -o errexit
    source /cluster/bin/jobsetup
    module purge
    module add ${HOME}/modules/fenics-2018-07-16

    # Add Ocellaris to Python's search path
    odir=${HOME}/src/Ocellaris
    export PYTHONPATH=${odir}:$PYTHONPATH

    # Run Ocellaris
    python3 ${odir}/scripts/orun.py ocellaris.inp --silent --timeout 1600

This requires some output settings in the Ocellaris input files:

.. code-block:: yaml

    output:
        # ... normal output settings

        # Restart friendly output settings
        stdout_enabled: yes
        flush_interval: 60 # (seconds, default 5)
        hdf5_write_interval: 50 # adjust to write every ten minutes or so
        hdf5_only_store_latest: yes  # to save disk space

The ``orun.py`` script has more options, but when running under SLURM some of
them are not needed as they are picked up from the environment.

.. code-block:: console

    $ python3 orun.py --help
    usage: orun [-h] [--ncpus NCPU] [--interval INTERVAL] [--pystuck]
                [--timeout TIMEOUT] [--restarts RESTARTS] [--silent]
                [--mpirun MPIRUN]
                input_file

    Start an Ocellaris simulation with a "babysitter" that watches the stdout
    stream and kills the simulation if no output is produced over an extended
    period of time (default 10 minutes / 600 seconds). If the simulation writes
    restart files at regular intervals then the babysitter can be made to restart
    the simulation from the latest restart file a number of times (default 2
    restarts of the same file). The reason for this babysitter is that there are
    difficult to debug problems (probably in PETSc) that causes the simulation to
    be stuck at 100% CPU utilisation. No backtrace is available on Ctrl+C / SIGINT
    which would be the case if there was a infinite loop in the Python code, so
    most likely the error exists in a C extension.

    positional arguments:
    input_file            Name of inputfile on YAML format

    optional arguments:
    -h, --help            show this help message and exit
    --ncpus NCPU, -n NCPU
                            Number of MPI processes. Not used when running under
                            SLURM where SLURM_NTASKS is used instead (default: 1)
    --interval INTERVAL, -i INTERVAL
                            Output interval in seconds (default: 10)
    --pystuck             Enable pystuck on the root MPI rank. Most likely will
                            not work since most hangs happen in C++ code (default:
                            False)
    --timeout TIMEOUT, -t TIMEOUT
                            Output timeout in seconds. After this period of
                            inactivity the simulation is killed (default: 600)
    --restarts RESTARTS, -r RESTARTS
                            Number of restarts of the same file (input or restart
                            file). Every time the simulation writes a new
                            savepoint the counter is reset (default: 2)
    --silent, -s          Do not relay stdout from Ocellaris (default: False)
    --mpirun MPIRUN       The mpirun executable (default: mpirun)

You could argue that finding the root cause of any PETSc MPI hangs would be
better than this hack to work around the problem, but I do not have time to
debug problems that happens after 50 hours of running on 48 CPUs somewhere
deep inside PETSc when the same routine has been called with more or less
similar matrices many thousand times before in the simulation without any
problems. It would be different if the hang was consistent and happened
earlier ... SORRY!


merge_xdmf_timeseries.py - join multiple XDMF files into a single time history
------------------------------------------------------------------------------

When postprocessing a simulation that has been restarted it can be 
inconvenient in programs such as Paraview that the time steps are spread
out over a number of XDMF files. This script merges such XDMF files into
one XDMF (and one HDF5) file that contains all the time steps. If restart
overlap the latest version of a timestep is written since this will be the
one that was used in the further simulation

Example:

.. code-block:: console

  $ rm merged.*
  $ python3 merge_xdmf_timeseries.py mysim.xdmf mysim_restarted_*.xdmf merged.xdmf

This will produce ``merged.xdmf`` and ``merged.h5``.

The script will not work for aribtrary XDMF files! It probably only works on
XDMF files produced by Ocellaris (and probably FEniCS DOLFIN with the same
XDMF configuration settings).

Others
------

These are not used much by me (Tormod Landet) and may hence have bitrotted and
could need some work to function as intended. Think of them more as examples to
start from if you need something similar, and not finished solutions.

- ``plot_reports.py`` - plot Ocellaris time step reports with matplotlib,
  optionally save a HTML report with the plots embedded in the file.

- ``plot_memory_usage.py`` - plot the memory usage for an Ocellaris simulation
  based on log file data. You must have specified ``output/show_memory_usage:
  yes`` in the input file to have the MAX RSS memory information available

- ``restart2vtk.py`` - take one function from an Ocellaris restart h5 file
  and export it as a true DG field to a ``*.vtk`` file. Currently only
  implemented for scalar DG2 fields, should be easy to extend to other element
  types. The binary VTK file writer may be buggy, the ASCII writer works.

- ``slice_to_numpy.py`` - read a restart file and extract a 2D slice of
  velocities and pressures which is stored as a numpy array on disk. Assumes
  that the simulation is 2D (ignores the z direction)

- *Various plotting scripts* - some use the newer ``ocellaris_post`` result
  file readers, some are from earlier times and implement result file parsers
  themselves (and should be updated).