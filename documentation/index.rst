About Ocellaris
###############

.. title:: Ocellaris: a mass-conserving higher-order DG FEM solver for free-surface flows

Ocellaris is a mass-conserving DG FEM solver for sharp-interface multiphase free-surface flows. Ocellaris can simulate water entry and exit of objects in ocean waves with accurate capturing of the force on the object and the behaviour of the free surface. Some examples of what Ocellaris can do, including videos of the results, are shown in the `Ocellaris Blog`_. Ocellaris is named after the `Amphiprion Ocellaris <https://en.wikipedia.org/wiki/Ocellaris_clownfish>`_ clownfish and is released under the Apache License, version 2.0.

.. figure:: https://www.ocellaris.org/figures/cylinder_in_waves.jpg
    :width: 50%
    :align: center
    :alt: Picture a cylinder in waves from an Ocellaris simulation

    Cylinder in waves; from Ocellaris via Paraview_ and Blender_.

Ocellaris is implemented in Python and C++ with FEniCS_ as the backend for the mesh and finite element assembly. PETSc_ is used for solving the resulting linear systems. The code is developed in Python and C++ on `Bitbucket <https://bitbucket.org/ocellarisproject/ocellaris>`_ by use of the Git version control system. If you want to contribute to Ocellaris, please read `the guide to contributing <https://www.ocellaris.org/programmers_guide/guidelines.html>`_. Ocellaris is automatically tested on `CircleCI <https://circleci.com/bb/ocellarisproject/ocellaris/tree/master>`_ and the current CI build status is |circleci_status|.


.. _sec_documentation_and_user_guide:

Documentation and user guide
############################

For help installing/running Ocellaris, please use the `issue tracker <https://bitbucket.org/ocellarisproject/ocellaris/issues>`_ and select ``Component = User help``. Please read this user guide first. We may start a user forum in the future if there are sufficient requests (let us know). Also, please let us know if there documentation is unclear (or wrong) in certain areas, or if you would like to see more documentation on a specific topic. This documentation will not replace a basic course in PDEs or CFD, but if you have taken those courses we hope you will be able to use Ocellaris for your purposes without any prior exposure to DG FEM.

.. toctree::
   :maxdepth: 3

   Ocellaris release notes and project news blog <blog/index>
   community_guidelines
   user_guide/user_guide
   programmers_guide/programmers_guide
   license
   contributors
   zreferences


.. _Ocellaris Blog: https://www.ocellaris.org/blog/
.. _FEniCS: https://fenicsproject.org/
.. _PETSc: https://www.mcs.anl.gov/petsc/
.. _Paraview: https://www.paraview.org/
.. _Blender: https://www.blender.org/
.. |circleci_status| image:: https://circleci.com/bb/ocellarisproject/ocellaris.svg?style=svg
    :target: https://circleci.com/bb/ocellarisproject/ocellaris

This documentation is written by Tormod Landet and the :doc:`Ocellaris contributors <contributors>` and is under the same :doc:`license <license>` as the rest of the Ocellaris source code repository contents. See :ref:`sec_citing` for how to reference Ocellaris in your research.
