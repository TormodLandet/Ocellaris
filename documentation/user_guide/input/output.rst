.. _inp_output:

Output control
==============

All the following parameters have sensible defaults and can be left out. The
output prefix can be useful to control in which directory the output files end
up. The final file name of all output files will be ``output_prefix +
file name``.

.. code-block:: yaml

    output:
        prefix: lid_driven_cavity_flow
        dolfin_log_level: warning
        ocellaris_log_level: info


.. describe:: stdout_enabled

    Log to standard out (console) in addition to the log file. Default on.

.. describe:: flush_interval

    Flush the log and stdout at regular intervals (given in seconds). Defaults
    to 0. A value of a couple of seconds will not impact performance a lot and
    will make it easier to follow simulations running on clusters which
    typically buffers output for quite long befor you see anything

.. describe:: log_on_all_ranks

    Write log files on all MPI ranks and not only on rank 0. This can be very
    helpfull when debugging a hanging simulation

.. describe:: stdout_on_all_ranks

    Write to standard out on all ranks. This is not useful for large MPI runs,
    but can be somewhat useful for small single machine runs. Default off.

.. describe:: log_enabled

    Write a log file. It is strongly recommended to leave this turned on, Many
    warnings can be written to the log and on the console they may scroll past
    too fast.

.. describe:: log_append_to_existing_file

    If a log file is found, append to this, do not remove it and start with an
    empty file. Defaults to on so that restarts preserve the log.

.. describe:: show_memory_usage

    This makes the log look ugly, but can help in discovering where memory is
    being allocated in case you want to trim some memory off your simulation.
    Default off.

.. describe:: hdf5_write_interval

    Write restart file every N time steps (default 0, never write restart file)

.. describe:: hdf5_only_store_latest

    Remove the previous restart file when the new one is finished writing to
    disk. Restart files can be large and you may only need the latest file.

.. describe:: save_restart_file_at_end

    Defaults to on, write a restart file when the simulation ends

.. describe:: xdmf_write_interval

    Write XDMF 3D plot files for Paraview and similar programs every N time
    steps. Defaults to 0, never write plot files. It is recommended to write
    XDMF to be able to visualize what is happening

.. describe:: xdmf_flush

    Flush the XDMF file after each write so that it can be opened while the
    simulation is running.

.. describe:: vtk_write_interval

    Write vtk files every N time steps. This writer is slow, but can handle
    some higher order fields (DG2) without interpolation to CG1 which is done
    in the XDMF writer currently. Not recommended unless you absolutely need
    this.

.. describe:: vtk_binary_format

    Defaults to off, the binary writer currently has a bug so use the ASCII
    writer for now or fix the bug.

.. describe:: solution_properties

    Compute and print properties such as divergence, courant number etc. This
    takes almost no time and is highly recommended

.. describe:: Co_lim

    Stop the simulation if the Courant number exceeds this value, default 1000.

.. describe:: plot_mesh

    Write the mesh to a separate plot file in the start of the simulation

.. describe:: plot_facet_regions

    Write the mesh with each facet region in a different colour to a separate
    plot file in the start of the simulation. Useful for checking boundary
    conditions.

.. describe:: plot_bcs

    Write the mesh with each boundary region in a different colour to a
    separate plot file in the start of the simulation. Useful for checking
    boundary conditions.

    One boundary region can consist of multiple facet regions. Facet regions
    are also not needed at all, the inside_code can be used to specify boundary
    regions, see :ref:`inp_boundary_conditions`.

