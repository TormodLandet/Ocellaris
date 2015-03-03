Ocellaris
=========

Ocellaris is a work in progress to make a continuous and discontinuous Galerkin FEM solver for 
multiphase free surface flows. The current goal of the project is to simulate water entry and 
exit of objects in ocean waves with accurate capturing of the force on the object and the 
behaviour of the free surface.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for numerics, mesh and 
finte element calculations.

.. _FEniCS: http://fenicsproject.org/

Ocellaris is named after the "Amphiprion Ocellaris" Clownfish and is written as part of a PhD
project at the University of Oslo.

Installation and running
------------------------

Ocellaris requires a full installation of FEniCS with the PETSc linear algebra backend. There is no
installation other than downloading the code and running::

  python -m ocellaris INPUTFILE.INP
  
with both the ``ocellaris`` Python package and the FEniCS and SciPy packages in the Python PATH 
(``dolfin``/``numpy``/``matplotlib`` etc). Ocellaris currently supports Python 2 only, not Python 3. 

To test the code there are some demos in the ``demos/`` directory. A complete input file is provided
for the well known lid driven cavity flow test case.

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015. Ocellaris is licensed under the Apache 2.0 license.
