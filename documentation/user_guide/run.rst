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
inspector is usefull when Ocellaris is running (for plotting log files) and
after finishing (plotting restart h5 files). All numbers printed on the screen
and in the log file when Ocellaris is running should be accessible in 
:ref:`the Ocellaris Inspector program <OcellarisInspector>`.

.. _meshio: https://github.com/nschloe/meshio
.. _gmsh: http://gmsh.info
.. _XDMF: http://www.xdmf.org
.. _Paraview: https://www.paraview.org
.. _HDFView: https://www.hdfgroup.org/downloads/hdfview/
.. _HDF5: https://www.hdfgroup.org
