Ocellaris containers
====================

This directory contains descriptions for how to build Docker and Singularity
containers for running Ocellaris. Using containers is by far the easiest way to
set up a working FEniCS environment with the necessary dependencies for running
Ocellaris


Docker
------

The easiest way to get started is to use the Docker container that is used for
running all the automated tests of Ocellaris. Currently this means running 
something like this::

    docker run -it trlandet/fenics-dev:py3_CI

Refer to the `documentation of Docker <https://docs.docker.com/>`_ for more
information. To install Ocellaris in the container, run something like::

    git clone https://bitbucket.org/trlandet/ocellaris.git
    cd ocellaris
    pip3 install --user .

You can also just add the Ocellaris source directory to the ``PYTHONPATH``
environmental variable. This can allow keeping the Ocellaris source outside the
container and is how Ocellaris is mainly developed (using Singularity instead of
Docker to run the container).


Singularity
-----------

Singularity, unlike Docker, does not require root access to start a container,
which means that it is available on shared computational resources where Docker,
which effectively gives the user root on the host machine, is not installed.

The supplied ``Singularity`` file adds Ocellaris to the Docker image and creates
a Singularity image from this combination. For more information, see the main '
Ocellaris documentation for installation advice and `the Singularity 
documentation <https://www.sylabs.io/docs/>`_.
