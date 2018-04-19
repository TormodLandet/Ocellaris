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

    docker run trlandet/fenics-dev:py3_CI

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
documentation <http://singularity.lbl.gov/>`_.


Building images
---------------

This description is meant for the developers of Ocellaris as a reminder since
this is not a task performed every day:

1. Log into a machine where you have root / can run Docker. At the Norwegian
   universities an `OpenStack IaaS service <http://www.uh-iaas.no>`_ can be
   used to create a virtual machine where you have root. When using this IaaS,
   each time you boot the virtual machine you need to mount the disk you have
   created for the storage of container image files, something like::
   
        sudo mount /dev/sdb /media/SingularityStorage/

2. Get the latest version of Ocellaris and go to the ``containers`` directory::

        # Either
        cd /media/SingularityStorage/images/Ocellaris/containers/
        git pull
        
        # Or
        cd /media/SingularityStorage/images
        git clone https://bitbucket.org/trlandet/ocellaris.git Ocellaris
        cd Ocellaris/containers/

3. Get the latest FEniCS dev-env Docker container::

        docker pull quay.io/fenicsproject/dev-env

4. Build FEniCS with some additional dependencies in Docker. This will take
   some time (between 10 minutes and 1 hour depending on hardware speed). Make
   sure there are no leftover large files (Singularity images etc) in the 
   directory since these will slow the process down due to being copied into
   the temporary build directory::

         # Run in the Ocellaris/containers directory
         docker build .

5. Figure out the IMAGE ID of the Docker image you just built::

        docker images

6. Tag the Docker image::

        docker tag XX-IMAGE-ID-XX trlandet/fenics-dev:py3_CI

7. Login and push the Docker image::

        docker login
        docker push trlandet/fenics-dev:py3_CI

8. Create a Singularity image (Singularity version 2.3, this changed in 2.4)::

        # Singularity version 2.3
        singularity create -s 3000 ocellaris.img
        sudo singularity bootstrap ocellaris.img Singularity
        
        # Singularity version 2.4
        sudo singularity build ocellaris.simg Singularity

9. Cleanup files (optional)::

        # List Docker containers
        docker ps -a
        # Delete Docker containers
        docker rm HASH1 HASH2 HASH3 ...
        
        # List Docker images
        docker images
        # Delete Docker images
        docker rmi HASH1 HASH2 HASH3 ...
        
        # Clean up Singularity's cache files
        sudo rm -r /root/.singularity
