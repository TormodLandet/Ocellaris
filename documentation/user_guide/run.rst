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
There exist a lot of programs to generate meshes that can be used with
Ocellaris. A Python program, `meshio-convert`_, can convert many formats
to FEniCS/DOLFIN/Ocellaris compatible XML format. A rather good free graphical
user interface for mesh generation is `gmsh`_. The latest version (>3) has
CAD capabilities and is used for several of the Ocellaris demos.

**2D/3D field postprocessing**:
Ocellaris can save results in XDMF_ format. There are several programs that
can postprocess such files. Paraview_ is a good option.


Custom postprocessor
~~~~~~~~~~~~~~~~~~~~

Ocellaris contains a custom postprocessor, the Ocellaris Inspector, which can
plot time step reports from ``*.log`` and restart files. The inspector can also
show the input and log files inside each restart file and 2D iso-line contours
like the location of the free surface in a 2D simulation. 

More advanced 2D and all 3D visualization is delegated to tools like Paraview_,
but for plotting the Courant number or total energy as a function of time step
for a finished simulation (log/restart file) or running simulation (log file)
the inspector can be quite handy. Using :kbd:`Control-R` to reload the results
gives near instant monitoring of a running Ocellaris simulation.

To run the Inspector you need to have a working installation of wxPython_
version 4.0 or later (tested with 4.0.0 beta 2). When this is installed you
can start the Ocellaris Inspector by running::

    python3 -m ocellaris_post.inspector myfile.log otherfile_endpoint_001.h5

This should show plots comparing the selected simulations. 

WxPython is not installed in Docker and Singularity containers for size
reasons, but the whole ``ocellaris_post`` Python package is written to not
depend on FEniCS or Ocellaris, so you can quite easily run it in a more
"standard" Python environment. The Ocellaris inspector currently supports
both Python 2 and Python 3 so you only need to install wxPython, numpy, h5py,
PyYAML and matplotlib in your favourite Python install (or venv) to use the 
Ocellaris Inspector. Most of these are probably installed allready in a
standard scientific Python installation.

For publication quality plots it is probaly best to use the ``ocellaris_post``
package from your own Python scripts instead of using Ocellaris Inspector.
Script examples can be found in the ``scripts/`` directory, though some of these
predate the ``ocellaris_post`` package. All results that are plotted in the
Inspector can be recreated by use of the :class:`ocellaris_post.Results` 
class:

.. autoclass:: ocellaris_post.Results

  .. attribute:: reports
        :annotation:
        
        A dictionary of report time series, the same as can be plotted in
        Ocellaris Inspector (Courant number, total mass, energy etc). 
        Corresponding x-axis data can be found in  :py:attr:`reports_x`.
    
  .. attribute:: reports_x
        :annotation:
        
        A dictionary of report time series x-axis data vectors. Same
        length as the corresonding entry in :py:attr:`reports`.
        
  .. attribute:: input
        :annotation:
        
        A plain dictionary conversion of the Ocellaris input file as
        returned by PyYAML.

.. _meshio-convert: https://github.com/nschloe/meshio
.. _gmsh: http://gmsh.info
.. _XDMF: http://www.xdmf.org
.. _Paraview: https://www.paraview.org
.. _HDFView: https://www.hdfgroup.org/downloads/hdfview/
.. _HDF5: https://www.hdfgroup.org
.. _wxPython: https://wxpython.org/
