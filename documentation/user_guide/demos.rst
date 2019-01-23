.. _demos:

Demos
=====

The demo input files can be found in the `Ocellaris repository
<https://bitbucket.org/trlandet/ocellaris/src/master/demos/>`_. Geometry files
(meshes) can be found in the ``datafiles`` subdirectory of the linked ``demos``
directory for the demos that use more complicated geometries.

The demos are automatically tested to check that they start properly, but since
some of the demos take quite some time to run, the demos are not run all the
way through. For this reason it may be that a change to Ocellaris can cause a
problem for one of the demos without anyone noticing. If you test a demo that
does not seem to work properly then please report an issue on the `Ocellaris
bug tracker <https://bitbucket.org/trlandet/ocellaris/issues>`_.


Selected demos
--------------


Flow around a clownfish
.......................

This demo shows how to create a simple Ocellaris 2D simulation. The gmsh
geometry file ``demos/datafiles/ocellaris.geo`` defines "physical regions" with
numeric IDs that can be referenced when defining boundary conditions in the
input file ``flow_around_ocellaris.inp``.

.. figure:: https://ocellarisproject.bitbucket.io/figures/flow_around_ocellaris.png
    :align: center
    :alt: Streamlines of the flow around a 2D clownfish

    The solution visualised in Paraview with an overlaid picture of the
    Ocellaris logo. The background color shows the distribution of the
    velocity magnitude, and the white lines show the stream lines of the
    converged steady state solution.


Dam break
.........

There are both 2D and 3D demos of the classic dam breaking two-phase flow
benchmark.

.. raw:: html

    <div class="figure align-center">
        <video controls loop autoplay>
            <source src="https://ocellarisproject.bitbucket.io/figures/dambreak2d.mp4" type="video/mp4">
            <source src="https://ocellarisproject.bitbucket.io/figures/dambreak2d.ogg" type="video/ogg">
            Your browser does not support the video tag.
        </video>
        <p class="caption">
            <span class="caption-text">
                Course mesh 2D dam break VOF simulation. Video from Paraview.
            </span>
        </p>
    </div>

The 2D simulation takes approximately 30 minutes to run on 1 CPU on a 2018
model laptop. Due to the course mesh it does not help using more CPUs. See
also :doc:`this blog entry <../blog/2018/11_green_water_movie>` for a 3D dam
breaking simulation that is a bit more impressive.


Taylor-Green
............

This is the standard Taylor-Green single phase 2D analytical solution. If you
want to run a convergence test then you should instead run the
``cases/convergence-taylor-green/convergence.py`` script which automatically
runs through a list of mesh resolutions and reports the convergence rates.

.. figure:: https://ocellarisproject.bitbucket.io/figures/taylor-green.png
    :align: center
    :alt: Taylor-Green velocity field

    The velocity field magnitude from the Taylor-Green demo

The simulation should only take a couple of seconds if the UFL has been
compiled to C++. On the first run this will happen, so this may take a bit
longer, but less than a minute.
