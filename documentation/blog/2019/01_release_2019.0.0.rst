New Ocellaris release, version 2019.0.0
=======================================

.. feed-entry::
    :author: Tormod Landet
    :date: 2019-01-10

A new version of Ocellaris has been released! This is a big milestone for the
project as all the essential features planed for Ocellaris are now implemented
and working 🎉 🎉 🙌

* Full 3D solver
* MPI parallel
* Two-phase flows
* Exactly mass conserving
* Higher order velocity approximation (quadratic polynomials in each cell)
* Stable and sharp air/water density transitions (factor 1000 density jump)

The new release has been tested successfully on several free surface flow
benchmarks on up to 64 CPUs on the UiO Abel HPC cluster and it gives very
satisfactory results. A paper on these tests is on the way. And, as always, a
set of shorter MMS tests such as Taylor-Green are run automatically on each
change in the code, so the basic flow solver should be in good shape.

.. figure:: https://ocellarisproject.bitbucket.io/figures/cylinder_in_waves_500.png
    :align: center
    :alt: Waves passing a cylinder, Ocellaris + Paraview + Blender + Cycles

    Waves passing a cylinder.
    Made with Ocellaris, Paraview_, Blender_ and Cycles_

.. _Paraview: https://www.paraview.org/
.. _Blender: https://www.blender.org/
.. _Cycles: https://www.cycles-renderer.org/

Ocellaris version 2019.0.0 is built on top of the latest FEniCS_ release as of
January 2019, `version 2018.1.0.r3 <https://quay.io/repository/fenicsproject/stable?tab=tags>`_.


Improved documentation
----------------------

The biggest feature new feature is the improved documentation. All relevant
input file parameters are now documented. A news feed/blog has been added to
show release notes such as this and other relevant information such as examples
of Ocellaris simulations.

Some of the more obscure parameters, which you do not want to use unless you
are deep into the code of Ocellaris itself, are not included in the
documentation on purpose. These parameters can still be found, either in the
source code itself (of course) or in the `complete list of input parameters
<https://bitbucket.org/ocellarisproject/ocellaris/src/master/ocellaris/input_file_schema.yml>`_
that is used to warn the user about misspellings in the input file. If there is
a parameter that you would like to see documented, or if something is unclear
in the documentation of a parameter, then please `file an issue on the bug
tracker <https://bitbucket.org/ocellarisproject/ocellaris/issues>`_.


Testing the new release
-----------------------

The release is available on `PyPi <https://pypi.org/project/ocellaris/2019.0.0/>`_,
but you will probably want to install it through either Docker or Singularity
unless you allready have an up to date FEniCS_ installation. The Docker method
should work on Linux, Mac OS X and Windows as long as Docker_ is installed.

More information about Docker and Singularity containers can be found in the
installation section of the Ocellaris documentation. To start a Docker
container for a quick test you can run::

    docker run -it trlandet/fenics-dev:2018.1.0.r3

The containers do not include Ocellaris, but contains all dependencies. To
install the 2019.0.0 version of Ocellaris inside the container you can run::

    pip3 install ocellaris==2019.0.0 --user

.. _FEniCS: https://fenicsproject.org/
.. _Docker: https://docs.docker.com/get-started/
