Running Ocellaris
-----------------


.. contents:: Contents
    :local:


Running a simulation
....................

Ocellaris is typically run from the command line with the name of an input file
as the first argument:

.. code-block:: sh

    ocellaris taylor-green.inp
    
You can optionally override parameters given on the input file:

.. code-block:: sh

    ocellaris taylor-green.inp \
        --set-input time/dt=0.1 \
        --set-input 'solver/velocity_function_space="CG"'  

You can see a summary of the command line options by running:

.. code-block:: sh

    ocellaris --help
    
Ocellaris will normally create a log file with the information that is also 
shown on screen. This will contain a listing of the input after modification
by the ``--set-input`` command line flag so that you can be sure to know
exactly what you did run when you look back on an old simulation.


Running a simulation on multiple CPUS with MPI
..............................................

The procedure for running in parallel with MPI is different depending on
whether you are running on an HPC (high performance cluster/computer), or on a
single machine. On a cluster---where you can run Ocellaris in parallel across
many machines---the method of running will also depend on whether your cluster
has FEniCS or Singularity installed. One of these must be installed for
Ocellaris to work.

**Not on HPC**

If you want to run in parallel on a single machine, then you will need to
specify the number of parallel processes to the ``mpirun`` command, and it
will launch and handle the parallel simulation. This also works inside a
Singularity or Docker container. To run on 8 CPUs:

.. code-block:: sh

    mpirun -np 8 ocellaris INPUTFILE.inp

You can see the number of CPUs that are being used in the Ocellaris log file
to verify that it is working as intended. The log file will also show the
number of mesh cells and function space unknowns per CPU.

**FEniCS installed on HPC**

If you have an installation of FEniCS at your local HPC cluster, and this
version of FEniCS is compiled for your clusters favourite MPI implementation,
then you can run Ocellaris just like any other MPI software (with the correct
``mpirun`` for the selected MPI implementation):

.. code-block:: sh

    # In a job-script for SLURM, Grid Engine, etc.
    mpirun ocellaris INPUTFILE.inp

The ``mpirun`` command (may also be called other names, check with your HPC
documentation) will typically figure out how many CPUs you have been allocated,
you normally do not have to specify this in the arguments to ``mpirun``. You
must typically write the above command inside a job-script and ask the HPC
scheduler to run the script for you with a given number of CPUs. How you write
the job-script and how you give it to the scheduler depends on your local HPC
machine, and should be documented by the HPC administrators.

You can also use :ref:`the orun.py script <script_orun>` to babysit your
simulation and make sure it is restarted if the cluster has a problem and the
simulation freezes (this can unfortunately happen with complex parallel file
system/MPI/Infiniband problems that are typically outside your control). You
replace the call to ``mpirun`` with ``python3 orun.py ...``. The ``orun``
script currently only supports SLURM, but it is relatively easy to extend it to
other schedulers. Send us a patch if you do!

**Singularity installed on HPC**

If you want to run on an HPC system where you do not have FEniCS available,
but you do have Singularity, then you can use an `MPICH ABI compatible
<https://www.mpich.org/abi/>`_ mpirun from outside Singularity to launch
your job:

.. code-block:: sh

    # With a scheduler
    mpirun singularity run ocellaris.img INPUTFILE.inp

    # Without a scheduler
    mpirun -np 16 singularity run ocellaris.img INPUTFILE.inp

You will need to use an MPICH ABI compatible ``mpirun``, such as Intel MPI or
MVAPICH2. You will probably also need to modify the Singularity image file so
that you include the appropriate Infiniband drivers. As an example, see
`this description <https://www.nsc.liu.se/support/singularity/mpi/>`_ or 
`this email <https://groups.google.com/a/lbl.gov/d/msg/singularity/50AhKYYQZVc/TAkx3oVhAQAJ>`_
for inspiration. By using a compatible MPI inside and outside the container it
should be just as fast as a local install (from other people's tests). This
has not been tested on Ocellaris simulations, let us know if you have some data
on this! Ocellaris is built with the standard MPICH library since this is what
is used by the FEniCS Docker images that we use as a basis.

Talk to your local HPC admin to figure out which libraries are necessary to
install inside the Singularity container. They may already have documents on
using Singularity if they have installed in on the HPC system. We cannot
document the exact steps to get it running here as HPC documentation is
unfortunately still very system specific. Hopefully tools like Singularity and
the push for reproducible science will improve this situation over time.


Restart files
.............

Ocellaris will by default save a restart file at the end of each simulation,
named something like ``SIMNAME_endpoint_00000XYZ.h5``. You can also configure
Ocellaris to write restart files at given intervals or supply a user code that
writes a restart file when given criteria are met. The restart file contains
the input file that was used along with a very detailed log and all active 
fields (velocity, pressure, density etc).

If you need to restart from the end of a simulation, for example to run a bit
further in time in case you set ``tmax`` a bit too short you can easily do this
by::

    ocellaris RESTART_FILE.h5 --set-input time/tmax=30.0 

You will probably want to use ``--set-input`` since it is inconvenient (but
certainly doable if you *really* want) to change the input description inside
the restart file.

If you want you can inspect the contents of a restart file, which is stored on
HDF5 format, by use of the graphical program HDFView_, or command line 
applications like ``h5ls`` and friends, see HDF5_ for more info. This is a
relatively easy way to find out the exact input file that was used, read the 
detailed log etc. There are no known programs that can plot the fields inside a
FEniCS DOLFIN HDF5 restart file (which is what Ocellaris uses).


Graphical user interface
........................

**Preprocessing:**
There are many programs that can generate simplical 2D or 3D meshes that are
compatible with Ocellaris. All mesh formats supported by `meshio`_, can be read
by Ocellaris. A rather good free graphical user interface for mesh generation is
`gmsh`_. The latest version (>3) has CAD capabilities and is used for several of
the Ocellaris demos.

**Postprocessing**:
Ocellaris can save results in XDMF_ format. There are several programs that
can postprocess such files. Paraview_ is one good option.
A custom post-processor, :ref:`OcellarisInspector`, also exist. It can be used
to plot residuals and other time series produced by the simulator. The 
inspector is useful when Ocellaris is running (for plotting log files) and
after finishing (plotting restart h5 files). All numbers printed on the screen
and in the log file when Ocellaris is running should be accessible in 
:ref:`the Ocellaris Inspector program <OcellarisInspector>`.

.. _meshio: https://github.com/nschloe/meshio
.. _gmsh: http://gmsh.info
.. _XDMF: http://www.xdmf.org
.. _Paraview: https://www.paraview.org
.. _HDFView: https://www.hdfgroup.org/downloads/hdfview/
.. _HDF5: https://www.hdfgroup.org


Controlling a running simulation
................................

If you run Ocellaris directly from a command line prompt on one CPU (very
unlikely) then you can type in the below commands directly. This is the reason
for the short command lines, to enable quick and dirty use during debugging of
the Ocellaris code on small toy examples.

If you are running Ocellaris in a cluster queue system and/or using MPI and
hence have no access to the interactive console of the root process you can give
commands in a command file. If your output prefix is such that the log file is
called ``mysim.log`` then the name of the command file is ``mysim.COMMANDS``. 
This file will not be created, you must make a new one every time. Ocellaris
will read the file to check for commands and then DELETE the file to avoid
running commands multiple times. The contained commands are processed at the end
of the current time step. Each command should be on a separate line in the file.

.. describe:: d

    Start a debug console - only use this when running interactively on 1 CPU!

.. describe:: f

    Flush open files to disk (this is done periodically, but can be forced)

.. describe:: p

    Plot field variables (legacy, avoid this).

.. describe:: r

    Write restart file

.. describe:: s

    Stop the simulation, changes input value ``time/tmax`` to the current time.

.. describe:: t

    Show timings (shows the table normally shown at the end of a simulation)

.. describe:: i a/b/c = 12.3

    Change the input file variable ``a/b/b`` to the value ``12.3``. The value
    will be evaluated as a Python expression.

.. describe:: w FORMAT

    Write current simulation state to 3D visuamization output file of the given
    FORMAT (must be one of ``vtk`` or ``xdmf``).

.. describe:: prof N

    Run the Python profiler for the next N time steps and then print the
    resulting profiling information. This is good for figuring out which routine
    is tanking longer than expected if the timing information from the ``t``
    command is not sufficient to understand the problem. Normally the PETSc
    Krylov solvers should take the majority of the time, but if not some 
    profiling may be necessary. Note: this profiles the Python code only, the
    C++ code will not show any details on which part is taking time.

**Example**: the following COMMAND file will save a restart file, plot fields to
XDMF and then stop the simulation::

    r
    w xdmf
    s

You can run something like this to easily create the COMMAND file::

    echo r > mysim.COMMANDS
    echo "w xdmf" >> mysim.COMMANDS
    echo s >> mysim.COMMANDS
