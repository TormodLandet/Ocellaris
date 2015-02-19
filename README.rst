Ocellaris
=========

A discontinuous Galerkin FEM solver for multiphase free surface flows. The current goal of the 
project is to simulate water entry and exit of objects in ocean waves.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for numerics, mesh and 
finte element calculations.

.. _FEniCS: http://fenicsproject.org/

Ocellaris is named after the "Amphiprion Ocellaris" Clownfish.

Installation and running
------------------------

Ocellaris requires a full installation of FEniCS with the PETSc linear algebra backend. There is no
installation other than downloading the code and running::

  python -m ocellaris INPUTFILE.INP
  
with both the ocellaris Python package and the FEniCS (dolfin/numpy etc) packages in the Python PATH.

To test the code there are some demos in the ``demos/`` directory. 

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015. Ocellaris is licensed under the Apache 2.0 license.
