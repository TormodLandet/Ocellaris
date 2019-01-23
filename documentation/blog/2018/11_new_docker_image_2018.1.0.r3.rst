New Docker image with FEniCS 2018.1.0.r3
========================================

.. feed-entry::
    :author: Tormod Landet
    :date: 2018-11-16

The Docker (and Singularity) container configurations have been updated
with the latest 2018.1.0.r3 release of FEniCS. The CircleCI test system is
now using this latest image and it is also available from `Docker Hub
<https://hub.docker.com/r/trlandet/fenics-dev/tags/>`_

.. cut::

More information can be found in the Ocellaris container documentation
(see :ref:`label-containers`). To start a container for a quick test you
can run::

    docker run -it trlandet/fenics-dev:2018.1.0.r3

The containers do not include Ocellaris, but contains all dependencies. To
install the latest version of Ocellaris inside the container you can run::

    pip3 install git+https://bitbucket.org/ocellarisproject/ocellaris.git --user
