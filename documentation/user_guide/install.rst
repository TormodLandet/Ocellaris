Installing Ocellaris
====================

Ocellaris is a Python package and it contains no modules that must be compiled before running. Some internal modules will be automatically compiled on the first program start-up by use of the FEniCS DOLFIN JIT compiler. This can take some time, up to 30 seconds. Subsequent runs of Ocellaris will use the pre-compiled modules.

The easiest way to install Ocellaris is by use of a Docker_ or Singularity_ container. Ocellaris is automatically CI_ tested after each change inside a Docker container, but the main development takes place using Singularity, and this is the recommended way to install Ocellaris.

.. contents:: Contents
    :local:


Dependencies
------------

**FEniCS**

Ocellaris depends on the FEniCS_ finite element framework---version 2018.1 or newer---with PETSc_ enabled (it is enabled by default). The most common ways ways to install FEniCS for use in Ocellaris are:

1) Use the provided Singularity container where FEniCS is pre-installed (**recommended**).

2) Use the provided Docker container where FEniCS is pre-installed.

3) Install from Debian or Ubuntu repositories, the following commands have been tested on Debian Buster:

  .. code-block:: bash

    apt-get install python3-lxml fenics python3-h5py python3-pip git

4) See the `FEniCS installation page`_ for details on other ways to install FEniCS.

**Python packages**

When you have installed FEniCS, the rest of the dependencies can be handled by the Python installation program, ``pip``. Simply running

.. code-block:: bash

  python3 -m pip install ocellaris
  
will install Ocellaris along with all dependencies (see :ref:`label-pip` for more details on using pip). For completeness, the Python package dependencies are also listed below:

* PyYAML
* h5py
* numpy
* matplotlib
* meshio >= 2.0.0
* raschii >= 1.0.2
* yschema >= 1.0.2

If you want to run the :ref:`postprocessing GUI <OcellarisInspector>`, you will additionally need to `install wxPython <https://wxpython.org/pages/downloads/>`_.


.. _label-containers:

Installation using containers
-----------------------------

**Singularity**

Ocellaris is developed (mostly) using Singularity_ containers. Singularity containers do not require root access to run, and are supported on (some) HPC clusters. You can either pull a pre-built image from SingularityHub,

.. code-block: bash

  singularity pull library://trlandet/default/ocellaris:2019.0.2

 which will leave you with an ``ocellaris.sif`` file in the current director, or you you can use the build description described in the ``Singularity`` file inside the  containers_ subdirectory of the Ocellaris source code. To build a Singularity image from the master branch of Ocellaris run:

.. code:: shell

  # Checkout the Ocellaris source code
  mkdir path/to/Ocellaris_source
  cd path/to/Ocellaris_source
  git clone https://bitbucket.org/ocellarisproject/ocellaris.git .
  # alternative: git clone git@bitbucket.org:ocellarisproject/ocellaris.git .

  # Build the Singularity image
  cd containers
  sudo singularity build ocellaris.img Singularity

This was tested with Singularity version 3.0.0. After building the image you can run Ocellaris from inside the newly created Singularity container::

  singularity run ocellaris.img INPUTFILE.INP

You can also get a shell session inside the image for testing and development::

  singularity shell ocellaris.img

**Docker**

Ocellaris is automatically CI_ tested after each change of the code by using the Docker container described in the containers_ subdirectory of the Ocellaris source code. You can see what this image is called on DockerHub in the CI config file `config.yml`_. Currently it is called ``trlandet/fenics-dev:py3_CI``, and you can launch your own Docker session using this container by running::

  docker run -it trlandet/fenics-dev:py3_CI

This will first download the image from Docker Hub and then show the standard FEniCS welcome message, but there will be some Ocellaris specific dependencies available in addition to a standard FEniCS distribution. The container does not include Ocellaris itself, but that is easily installed with pip now that all the dependencies are present::

  python3 -m pip install --user ocellaris


.. _label-pip:

Installation using pip
----------------------

You can install the latest stable version of Ocellaris by running::

    python3 -m pip install ocellaris

To install the master version you must check out the git repository and install from there:

.. code:: shell

  git clone https://bitbucket.org/ocellarisproject/ocellaris.git
  # alternative: git clone git@bitbucket.org:ocellarisproject/ocellaris.git
  cd ocellaris
  python3 -m pip install .

If Ocellaris is installed via ``pip``, then the ``ocellaris`` command will be available, otherwise you can add the source directory of Ocellaris to the Python module search path manually and add an alias::

    alias ocellaris="python3 -m ocellaris"

Ocellaris depends on a working installation of FEniCS, compiled with support for PETSc, and some additional Python packages like PyYAML and h5py. Ocellaris will inform you about any missing packages when you run it for the first time. FEniCS is not currently pip-installable (as of February 2019), so it can be slightly hard to install all prerequisites. The recommended way, which should always work, is to use the same installation as on the automated test system—running in a container. You can also install from ``deb`` packages or conda, see the `FEniCS web pages <https://fenicsproject.org/download/>`_. FEniCS power users can compile and install on their own, but the main author of Ocellaris once used a week to get a fast and working FEniCS `installed on an old-ish HPC cluster <https://bitbucket.org/trlandet/fenics-on-abel>`_, it is not always easy ...


Verifying the installation
--------------------------

The automated testing configuration file, `config.yml`_, shows you how to run the Ocellaris test on your own machine to verify that the installation works as intended (the test commands are also :ref:`shown here <label-running-tests>`).

You can also run one of the demos:

.. code:: shell

  cd path/to/Ocellaris_source
  cd demos
  # Use one of these commands to run an Ocellaris simulation
  ocellaris taylor-green.inp
  python3 -m ocellaris taylor-green.inp


.. _FEniCS: https://www.fenicsproject.org/
.. _PETSc: https://www.mcs.anl.gov/petsc/
.. _`FEniCS installation page`: https://fenics.readthedocs.io/en/latest/installation.html
.. _config.yml: https://bitbucket.org/ocellarisproject/ocellaris/src/master/.circleci/config.yml
.. _containers: https://bitbucket.org/ocellarisproject/ocellaris/src/master/containers
.. _Singularity: https://www.sylabs.io/singularity/
.. _Docker: https://www.docker.com/
.. _CI: https://circleci.com/bb/ocellarisproject/ocellaris/tree/master
