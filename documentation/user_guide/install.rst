Installing Ocellaris
--------------------

Ocellaris is a Python package and it contains no modules that must be compiled
before running. Some internal modules will be compiled on the first program
startup by use of the FEniCS DOLFIN JIT compiler. This can take some time.
Subsequent runs of Ocellaris will use the precompiled modules.


.. contents:: Contents
    :local:


.. _label-containers:

Installation using containers
.............................

The easiest way to install Ocellaris is by use of a Docker_ or Singularity_
container. Ocellaris is CI_ tested using the Docker container described in
the ``containers/`` subdirectory of the Ocellaris source code. The test
procedures (Linux shell commands) describe the exact commands used to
install and run Ocellaris tests and they are a good place to start, see the
`config.yml`_ file in the ``.circleci/`` subdirectory of the Ocellaris source
code for the details.  

Ocellaris is developed (mostly) using Singularity_ containers. You can either
convert the Docker CI container or use the one described in the ``Singularity``
file inside the ``containers/`` subdirectory of the Ocellaris source code.
To create a Singularity image from the file "Singularity" run::

  cd path/to/Ocellaris_source
  cd containers
  singularity build ocellaris.img Singularity

You can now run Ocellaris from inside the newly created Singularity container::

  singularity run ocellaris.img INPUTFILE.INP

The Singularity image is based on the Docker image that is used by the Ocellaris
automated testing environment, see `config.yml`_ for up to date details about
which Docker image is used. You can use this Docker image to run Ocellaris as
well, but you will then have to install Ocellaris yourself inside the container
using the ``pip3`` command shown below.


Installation using pip
......................

Before running Ocellaris you must ensure that the ``ocellaris`` Python package
is on the Python search path. This is most easily done by running::

    pip3 install .
   
in the root directory of the source code. If the package is installed via
``pip`` then the ``ocellaris`` command will be available, otherwise you can
add the source directory to the Python module search path and add an alias::

    alias ocellaris="python3 -m ocellaris"

Ocellaris depends on an installation of FEniCS, compiled with support for
PETSc, and some additional Python packages like PyYAML and h5py. Ocellaris will
inform you about any missing packages when you run it for the first time.

The Ocellaris version will be available on PYPI for installation may lag 
significantly behind the master version. This is the version you get from
``pip install ocellaris``. To get the latest master version you can 
download the source code manually from  `the Ocellaris Bitbucket git repository
<https://bitbucket.org/trlandet/ocellaris/src>`_. You can get the source code
and install it by running::

  git clone https://bitbucket.org/trlandet/ocellaris.git
  cd ocellaris
  pip3 install .

FEniCS, which Ocellaris is built on top of, is not currently pip-installable
as of December 2018, so it can be slightly hard to install all prerequisites.
The recommended way which should always work is to use the same installation as
on the automated test system—running in a container—or using the same
installation procedure as used in the containers, see the container section
above for more info. 


.. _config.yml: https://bitbucket.org/trlandet/ocellaris/src/master/.circleci/config.yml
.. _Singularity: http://singularity.lbl.gov/
.. _Docker: https://www.docker.com/
.. _CI: https://circleci.com/bb/trlandet/ocellaris/tree/master


