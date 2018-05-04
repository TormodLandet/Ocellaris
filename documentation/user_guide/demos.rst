Demos
=====

The demo input files can be found in the `Ocellaris repository
<https://bitbucket.org/trlandet/ocellaris/src/master/demos/>`_. Geometry files
(meshes) can be found in the ``datafiles`` subdirectory of the linked ``demos``
directory.


Flow around a clownfish
-----------------------

This demo shows how to create a simple Ocellaris 2D simulation. The gmsh
geometry file ``demos/datafiles/ocellaris.geo`` defines "physical regions" with
numeric IDs that can be referenced when defining boundary conditions in the
input file ``flow_around_ocellaris.inp``.

.. figure:: https://trlandet.bitbucket.io/figures/flow_around_ocellaris.png
    :align: center
    :alt: Streamlines of the flow around a 2D clownfish
        
    The solution visualised in Paraview with an overlaid picture of the
    Ocellaris logo. The background color shows the distribution of the
    velocity magnitude, and the white lines show the stream lines of the
    converged steady state solution.

TODO: explain some details of the input file?


Dam break
---------

TODO: document this and include a figure


Taylor-Green
------------

TODO: document this and include a figure


Wave Tank
---------

TODO: document this and include a figure
