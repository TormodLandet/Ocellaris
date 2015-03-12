Ocellaris
=========

Ocellaris is a work in progress to make a continuous and discontinuous Galerkin FEM solver for 
multiphase free surface flows. The current goal of the project is to simulate water entry and 
exit of objects in ocean waves with accurate capturing of the force on the object and the 
behaviour of the free surface.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for numerics, mesh and 
finte element calculations.

.. _FEniCS: http://fenicsproject.org/

Ocellaris is named after the `Amphiprion Ocellaris <http://en.wikipedia.org/wiki/Ocellaris_clownfish>`_
clownfish and is written as part of a PhD project at the University of Oslo.

.. image:: http://trlandet.bitbucket.org/ocellaris/_static/ocellaris_mesh_521.png
    :align: center
    :alt: Picture of Ocellaris

Installation and running
------------------------

Ocellaris requires a full installation of FEniCS_ with the PETSc linear algebra backend. There is no
installation other than downloading the code and running the following command with both the Ocellaris
Python package and the FEniCS and SciPy packages in the Python PATH (dolfin/numpy/matplotlib etc)::

  python -m ocellaris INPUTFILE.INP
  
Ocellaris currently supports Python 2 only, not Python 3. 

To test the code there are some demos in the ``demos/`` directory. A complete input file is provided
for the well known lid driven cavity flow test case. More information can be found in the documentation.

Documentation
-------------

.. TOC_STARTS_HERE  - in the Sphinx documentation a table of contents will be inserted here 

The documentation can be found on the `Ocellaris web page <http://trlandet.bitbucket.org/ocellaris/>`_.

.. TOC_ENDS_HERE

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015. Ocellaris is licensed under the Apache 2.0 license.
