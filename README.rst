Ocellaris
=========

Ocellaris is a work in progress to make a mass conserving DG FEM solver for sharp interface
multiphase free surface flows. The current goal of the project is to simulate water entry and 
exit of objects in ocean waves with accurate capturing of the force on the object and the 
behaviour of the free surface.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for numerics, mesh and 
finite element calculations.

.. _FEniCS: https://fenicsproject.org/

Ocellaris is named after the `Amphiprion Ocellaris <https://en.wikipedia.org/wiki/Ocellaris_clownfish>`_
clownfish and is written as part of a PhD project at the University of Oslo.

.. figure:: https://trlandet.bitbucket.io/ocellaris/_static/ocellaris_mesh_521.png
    :align: center
    :alt: Picture of Ocellaris
    
    `About this image <https://trlandet.bitbucket.io/ocellaris/logo.html>`_

Installation and running
------------------------

Ocellaris requires a full installation of FEniCS_ with the PETSc linear algebra backend. You can
install the dependecies yourself (you need at least dolfin, h5py, matplotlib and PyYAML) and then
install ocellaris to somewhere in the python module search path. You can then run::

  python -m ocellaris INPUTFILE.INP

You can also install using the preliminary support for Singularity containers::

  # Bootstrap a Singularity image from the file "Singularity" which
  # is located in the root of the Ocellaris git repository
  singularity create -s 3000 ocellaris.img
  sudo singularity bootstrap ocellaris.img Singularity

  # Run bash inside the newly created Singularity container
  singularity run -H /some/empty/dir/to/use/as/home ocellaris.img

  # Or, just run Ocellaris directly, the container exits along with Ocellaris
  singularity run -H /some/empty/dir/to/use/as/home ocellaris.img -c "ocellaris INPUTFILE.INP"
  
To test the code there are some demo input files in the ``demos/`` directory. Complete input files along
with driver scripts are provided for several of the normal benchmark cases like Kovasznay flow and the
Taylor-Green vortex in the ``cases/`` directory. More information can be found in the documentation which
also contains an (incomplete) description of the input file format.

Please feel free to test Ocellaris, but please keep in mind:

- Ocellaris is in a state of constant development and does not have a stable API or input file format
- Ocellaris supports Python 2 only, not Python 3 (currently, hopefully this will change some day).
- This is a research project, do not expect anything to work properly without testing it thoroughly first!
- Documentation has not been a priority, sorry!

Documentation
-------------

.. TOC_STARTS_HERE  - in the Sphinx documentation a table of contents will be inserted here 

The documentation can be found on the `Ocellaris web page <https://trlandet.bitbucket.io/ocellaris/>`_.

.. TOC_ENDS_HERE

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015-2017. Ocellaris is licensed under the Apache 2.0 license.
