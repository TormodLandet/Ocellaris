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
