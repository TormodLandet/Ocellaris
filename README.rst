Ocellaris
=========

Ocellaris is a mass conserving DG FEM solver for sharp interface multiphase
free surface flows. Ocellaris can simulate water entry and exit of objects in
ocean waves with accurate capturing of the force on the object and the
behaviour of the free surface. Some examples of what Ocellaris can do,
including videos of the results, are shown in the `Ocellaris Blog`_.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for the
mesh and finite element assembly. PETSc_ is used for solving the resulting
linear systems.

.. contents:: Quick start

.. _Ocellaris Blog: https://trlandet.bitbucket.io/ocellaris/blog/
.. _FEniCS: https://fenicsproject.org/
.. _PETSc: https://www.mcs.anl.gov/petsc/

Ocellaris is named after the `Amphiprion Ocellaris <https://en.wikipedia.org/wiki/Ocellaris_clownfish>`_
clownfish and is written as part of a PhD project at the University of Oslo.

.. figure:: https://trlandet.bitbucket.io/figures/ocellaris_outlined_500.png
    :align: center
    :alt: Picture of an Ocellaris clownfish in a triangulated style

Installation and running
------------------------

Ocellaris requires a full installation of FEniCS_ with the PETSc linear algebra
backend. You can install the dependecies yourself (you need at least dolfin,
h5py, matplotlib and PyYAML), but the easiest way by far is to use a
preconfigured Singularity or Docker container. More information on these and
installation in general can be found in the `user guide`_.

When Ocellaris is installed you can run the solver with an Ocellaris input
file::

  ocellaris INPUTFILE.INP

Example input files can be found in the ``demos/`` sub-directory of the
Ocellaris source code and a description of the Ocellaris input file format and
the possible input parameters is given in the `user guide`_.

.. _user guide: https://trlandet.bitbucket.io/ocellaris/user_guide/user_guide.html

First steps
~~~~~~~~~~~

To test the code there are some demo input files in the ``demos/`` directory.
Complete input files along with driver scripts are provided for several of the
standard benchmark cases like Kovasznay flow and the Taylor-Green vortex in the
``cases/`` directory. More information can be found in the documentation which
also contains a description of the input file format.

Please feel free to test Ocellaris, but please keep in mind:

- Ocellaris is in a state of constant development
- Ocellaris is tested with FEniCS Version 2018.1. Earlier versions will NOT
  work, later version may possibly work.
- This is an ongoing research project, do not expect results to be correct
  without proper validation!


Documentation
-------------

.. TOC_STARTS_HERE  - in the Sphinx documentation a table of contents will be inserted here 

The documentation can be found on the `Ocellaris web page <https://trlandet.bitbucket.io/ocellaris/index.html#documentation>`_.

.. TOC_ENDS_HERE

Development
-----------

Ocellaris is developed in Python and C++ on `Bitbucket <https://bitbucket.org/trlandet/ocellaris>`_
by use of the Git version control system. If you are reading this on github,
please be aware that you are seeing a mirror that could potentially be months
out of date. The github mirror is only updated sporadically—to trigger new
Singularity and Docker Hub container builds. All pull requests and issues
should go to the Bitbucket repository.

Ocellaris is automatically tested on `CircleCI <https://circleci.com/bb/trlandet/ocellaris/tree/master>`_
and the current CI build status is |circleci_status|.

.. |circleci_status| image:: https://circleci.com/bb/trlandet/ocellaris.svg?style=svg&circle-token=886a679594f958395d69c0720b04c4d88056f49d
  :target: https://circleci.com/bb/trlandet/ocellaris/tree/master

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015-2019. Ocellaris is licensed under
the Apache 2.0 license, a permissive free software license compatible with
version 3 of the GNU GPL. See `License of Ocellaris`_ for the details.

.. _`License of Ocellaris`:  https://trlandet.bitbucket.io/ocellaris/license.html
