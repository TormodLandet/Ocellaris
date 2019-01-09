.. _inp_probes:

Probes
======

Probes can be used to ease postprocessing and to extract numbers directly from
the simulation instead of having to write XDMF files which are large and may
be somewhat interpolated from the "real" value inside the simulation. Often
XDMF and probes are used at the same time, probes write their values every time
step and XDMF may be written every 100 time steps for disk saving purposes.

The probes list contains dictionaries describing each probe. An example can be
seen below under :ref:`inp_probes_pointprobe`. Some common configuration
options for all the probes are:

.. describe:: name

    Name of the probe, often used in the name of the output file

.. describe:: type

    The probe type. One of

    * IsoSurface
    * LineProbe
    * PlaneProbe
    * PointProbe

.. describe:: enabled

    Default on, you can turn of a probe if you want

.. describe:: file_name

    Where to save the data. Normally a default value based on the name is used

.. describe:: write_interval

    How often to write the probe results, default 1, every time step

.. describe:: custom_hook

    When to run the probe. Normally run at the end oe each time step, but you
    can give the name of a hook here instead (typically
    ``MultiPhaseModelUpdated``). Not all probes support this.


IsoSurface
----------

Only implemented in 2D for now. Outputs the position of an ISO surface. Can be
used to save the free surface. Another way (for 3D) is to use XDMF output and
create a Contour in Paraview as a post-processing step.

.. describe:: field

    Name of the field to study, e.g., ``c``

.. describe:: value

    The value of the field at the ISO surface, e.g. ``0.5`` for the free
    surface in a VOF simulation


LineProbe
---------

.. describe:: field

    Name of the field to study, e.g., ``u0``

.. describe:: startpos

    A list of numbers, the coordinates of a point in the domain

.. describe:: endpos

    A list of numbers, the coordinates of a point in the domain

.. describe:: Npoints

    Number of probing points along the line segment from startpos to endpos


PlaneProbe
----------

Saves an XDMF plot file of the specified field intersected by a plane.
Sometimes it can be useful to have 2D slices of 3D simulations since the 2D
slices are smaller in size and can be written more often without too much IO.

.. describe:: field

    Name of the field to study, e.g., ``u0``. You can also give a list,
    ``[u0, u1, u2]``, but the functions in the list must share the same
    function space (most likely DG2 in this case)

.. describe:: plane_point

    A list of numbers, the coordinates of a point on the plane

.. describe:: plane_normal

    A list of numbers, the normal direction of the plane

.. describe:: xlim, ylim, zlim

    Lists of two numbers specifying limits to the extents of the plane. By
    default the plane is as large as the intersection with the 3D mesh allows.


.. _inp_probes_pointprobe:

PointProbe
----------

Probe one or more fields in given points

.. describe:: probe_points

    A list of function names and the coordinates of the points to probe. The
    name of each probe must also be given so that you can figure out which
    value belongs to which point. See example below for the syntax


.. code-block:: yaml

    probes:
        -   name: pressure_probes
            enabled: yes
            type: PointProbe
            probe_points:
            -   ['p', 'probe1', -1, 0.5, 0.3]
            -   py$ ['p', 'probe2', L - 2, 1e-3, 1e-3]
            -   py$ ['c', 'cprobe', L - 2, 1e-3, 1e-3]
