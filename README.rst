Ocellaris
=========

Ocellaris is a work in progress to make a mass conserving DG FEM solver for sharp interface
multiphase free surface flows. The current goal of the project is to simulate water entry and 
exit of objects in ocean waves with accurate capturing of the force on the object and the 
behaviour of the free surface.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for numerics, mesh and 
finite element calculations.

.. contents:: Quick start

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
install the dependecies yourself (you need at least dolfin, h5py, matplotlib and PyYAML) or use 
a preconfigured Singularity or Docker container, see below.

To install Ocellaris in an environment where the prerequisites are present run::

  cd path/to/ocellaris/source
  pip3 install .

You should now be ready to run the solver::

  ocellaris INPUTFILE.INP

If you did not pip install the ocellaris package, but just put it somewhere in the Python path,
you may need to instead run::

  python3 -m ocellaris INPUTFILE.INP

Singularity and Docker
~~~~~~~~~~~~~~~~~~~~~~

You can also install using the preliminary support for Singularity containers. First, create a Singularity image from the file "Singularity" which is located in the root of the Ocellaris git repository::

  cd path/to/ocellaris/source
  singularity build ocellaris.img Singularity

You can now run Ocellaris from inside the newly created Singularity container::

  singularity run ocellaris.img INPUTFILE.INP

The Singularity image is based on the Docker image that is used by the Ocallaris automated testing
environment, see `setup.yml <https://bitbucket.org/trlandet/ocellaris/src/master/.circleci/config.yml>`_
for up to date details about which Docker image is used. You can use this Docker image to run
Ocellaris as well, but you will then have to install Ocellaris yourself inside the container,
see details in the ``setup.yml`` or ``Singularity`` files, it is basically just ``pip install``.

First steps
~~~~~~~~~~~
  
To test the code there are some demo input files in the ``demos/`` directory. Complete input files along
with driver scripts are provided for several of the standard benchmark cases like Kovasznay flow and the
Taylor-Green vortex in the ``cases/`` directory. More information can be found in the documentation which
also contains a (currently incomplete) description of the input file format.

Please feel free to test Ocellaris, but please keep in mind:

- Ocellaris is in a state of constant development and does not have a stable API or input file format. 
- Ocellaris supports Python 3 only (after October 2017). There is an old Python2 branch ``legacy-py2``
  that may or may not work. Dolfin with pybind11 Python3 wrappers is required (master version, still
  unreleased as of 2017-10-17).
- This is a research project, do not expect anything to work properly without testing it thoroughly first!
- Documentation has not been a big priority, sorry! See the demos, regression tests etc in the code 
  repository in addition to the documentation and of course the source code itself.

Documentation
-------------

.. TOC_STARTS_HERE  - in the Sphinx documentation a table of contents will be inserted here 

The documentation can be found on the `Ocellaris web page <https://trlandet.bitbucket.io/ocellaris/>`_.

.. TOC_ENDS_HERE

Development
-----------

Ocellaris is developed in Python and C++ on `Bitbucket <https://bitbucket.org/trlandet/ocellaris>`_ by use
of the Git version control system. If you are reading this on github, please be aware that you are seeing a
mirror that could potentially be months out of date. The github mirror is only updated sporadically—to trigger
new Singularity Hub container builds. All pull requests and issues should go to the Bitbucket repository.

Ocellaris is automatically tested on `CircleCI <https://circleci.com/bb/trlandet/ocellaris/tree/master>`_  
and the current CI build status is |circleci_status|.

.. |circleci_status| image:: https://circleci.com/bb/trlandet/ocellaris.svg?style=svg&circle-token=886a679594f958395d69c0720b04c4d88056f49d
  :target: https://circleci.com/bb/trlandet/ocellaris/tree/master

Copyright and license
---------------------

Ocellaris is copyright Tormod Landet, 2015-2017. Ocellaris is licensed under the Apache 2.0 license, a 
permissive free software license compatible with version 3 of the GNU GPL.
