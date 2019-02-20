Installing Ocellaris
--------------------

Ocellaris is a Python package and it contains no modules that must be compiled before running. Some internal modules will be automatically compiled on the first program start-up by use of the FEniCS DOLFIN JIT compiler. This can take some time, up to 30 seconds. Subsequent runs of Ocellaris will use the pre-compiled modules.


.. contents:: Contents
    :local:


.. _label-containers:

Installation using containers
.............................

The easiest way to install Ocellaris is by use of a Docker_ or Singularity_ container. Ocellaris is automatically CI_ tested after each change of the code by using the Docker container described in the containers_ subdirectory of the Ocellaris source code. You can see what this image is called on DockerHub in the CI config file `config.yml`_. Currently it is called ``trlandet/fenics-dev:py3_CI``, and you can launch your own Docker session using this container by running::

  docker run -it trlandet/fenics-dev:py3_CI

This will first download the image from Docker Hub and then show the standard FEniCS welcome message, but there will be some Ocellaris specific dependencies available in addition to a standard FEniCS distribution. The container does not include Ocellaris itself, but that is easily installed with pip now that all the dependencies are present::

  pip3 install --user ocellaris

Ocellaris is developed (mostly) using Singularity_ containers. Singularity containers do not require root access to run and are supported on (some) HPC clusters. You can either convert the Docker CI container or use the Singularity build description described in the ``Singularity`` file inside the  containers_ subdirectory of the Ocellaris source code. To build a Singularity image from the master branch of Ocellaris run:

.. code:: shell

  # Checkout the Ocellaris source code
  mkdir path/to/Ocellaris_source
  cd path/to/Ocellaris_source
  git clone git@bitbucket.org:ocellarisproject/ocellaris.git .

  # Build the Singularity image
  cd containers
  singularity build ocellaris.img Singularity

This was tested with Singularity version 3.0.0. After building the image you can run Ocellaris from inside the newly created Singularity container::

  singularity run ocellaris.img INPUTFILE.INP

You can also get a shell session inside the image for testing and development::

  singularity shell ocellaris.img

The automated testing configuration file, `config.yml`_, shows you how to run the Ocellaris test on your own machine to verify that the installation works as intended. You can also run one of the demos:

.. code:: shell

  cd path/to/Ocellaris_source
  cd demos
  # Use one of these commands to run an Ocellaris simulation
  ocellaris taylor-green.inp
  python3 -m ocellaris taylor-green.inp

Installation using pip
......................

You can install the latest stable version of Ocellaris by running::

    pip3 install ocellaris

To install the master version you must check out the git repository and install from there:

.. code:: shell

  git clone https://bitbucket.org/ocellarisproject/ocellaris.git
  # alternative: git clone git@bitbucket.org:ocellarisproject/ocellaris.git
  cd ocellaris
  pip3 install .

If Ocellaris is installed via ``pip`` then the ``ocellaris`` command will be available, otherwise you can add the source directory of Ocellaris to the Python module search path manually and add an alias::

    alias ocellaris="python3 -m ocellaris"

Ocellaris depends on a working installation of FEniCS, compiled with support for PETSc, and some additional Python packages like PyYAML and h5py. Ocellaris will inform you about any missing packages when you run it for the first time. FEniCS is not currently pip-installable (as of February 2019), so it can be slightly hard to install all prerequisites. The recommended way, which should always work, is to use the same installation as on the automated test system—running in a container. You can also install from ``deb`` packages or conda, see the `FEniCS web pages <https://fenicsproject.org/download/>`_. FEniCS power users can compile and install on their own, but the main author of Ocellaris once used a week to get a fast and working FEniCS `installed on an old-ish HPC cluster <https://bitbucket.org/trlandet/fenics-on-abel>`_, it is not always easy ...

.. _config.yml: https://bitbucket.org/ocellarisproject/ocellaris/src/master/.circleci/config.yml
.. _containers: https://bitbucket.org/ocellarisproject/ocellaris/src/master/containers
.. _Singularity: https://www.sylabs.io/singularity/
.. _Docker: https://www.docker.com/
.. _CI: https://circleci.com/bb/ocellarisproject/ocellaris/tree/master


