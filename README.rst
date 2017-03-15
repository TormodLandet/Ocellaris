Ocellaris
=========

Ocellaris is a work in progress to make a continuous and discontinuous Galerkin FEM solver for 
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

Ocellaris requires a full installation of FEniCS_ with the PETSc linear algebra backend. There is no
installation other than downloading the code and running the following command with both the Ocellaris
Python package and the FEniCS and SciPy packages in the Python PATH (dolfin/numpy/matplotlib etc)::

  python -m ocellaris INPUTFILE.INP
  
To test the code there are some demos in the ``demos/`` directory. Complete input files are provided
for several of the normal benchmark cases like lid driven cavity flow and the Taylor-Green vortex. More
information can be found in the documentation which also contains a description of the input file format.

Please feel free to test Ocellaris, but please keep in mind:

- Ocellaris is in a state of constant development and does not have a stable API or input file format
- Ocellaris supports Python 2 only, not Python 3 (currently).
- This is a research project, do not expect anything to work properly without testing it thoroughly first!

Documentation
-------------

.. TOC_STARTS_HERE  - in the Sphinx documentation a table of contents will be inserted here 

The documentation can be found on the `Ocellaris web page <https://trlandet.bitbucket.io/ocellaris/>`_.

.. TOC_ENDS_HERE

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015-2017. Ocellaris is licensed under the Apache 2.0 license.