.. _OcellarisInspector:

Ocellaris Inspector
===================

Ocellaris contains a custom postprocessor, the Ocellaris Inspector, which can
plot time step reports from log and restart files. The inspector can also show
the input and log files inside each restart file and has a number of more 
specialized plotting functionality for 2D iso-lines and quasi-static 
analyses. 

More advanced 2D and all 3D visualization is delegated to tools like Paraview_,
but for plotting the Courant number or total energy as a function of time step
for a finished simulation (log/restart file) or running simulation (log file)
the inspector can be quite handy. Using :kbd:`Control+R` to reload the results
gives near instant monitoring of a running Ocellaris simulation.

.. contents:: Contents
    :local:

Getting started
---------------

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


.. figure:: https://ocellarisproject.bitbucket.io/figures/inspector_courant.png
    :align: center
    :alt: The Ocellaris Inspector GUI - plotting Co = Courant number
    
    The Ocellaris Inspector showing the Courant number of a running simulation 


Opening files
-------------

You can specify files on the command line, click the "Open" button in the
"Setup" tab, drag and drop files into the program or use :kbd:`Control+O` to
bring up a file opening dialog.

You can also load log files from simulations running on a HPC cluster system
by clicking the "Open running simulation" button in the "Setup" tab and giving
the host name and mount directories for the HPC cluster login node. You **must**
have enabled password-less SSH login to the head node, otherwise the cluster
connector will not work

.. figure:: https://ocellarisproject.bitbucket.io/figures/cluster_connector.png
    :align: center
    :alt: The cluster connection GUI

    The cluster connection GUI

You must also have mounted the cluster home directory somewhere on the local
machine; via ``sshfs`` or other means. When the connection is established you
can select which simulations you would like to open and press the "Load" button.

Only SLURM clusters are currently supported, but just a tiny bit of glue code is
necessary to implement support for other queue systems; please submit a patch or
get in touch if you need such support. Loading running simulations from the
current machine should also be easy and may be implemented in the future.


Scripting
---------

For publication quality plots it is probaly best to use the ``ocellaris_post``
package from your own Python scripts instead of using Ocellaris Inspector.
Script examples can be found in the ``scripts/`` directory, though some of
these predate the ``ocellaris_post`` package. All results that are plotted in
the Inspector can be recreated by use of the :class:`ocellaris_post.Results`
class. Open using the :func:`ocellaris_post.open_results` function:


.. autofunction:: ocellaris_post.open_results

.. autofunction:: ocellaris_post.read_yaml_input_file

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

.. _Paraview: https://www.paraview.org
.. _wxPython: https://wxpython.org
