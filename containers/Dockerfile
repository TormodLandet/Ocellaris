# Docker image used to test Ocellaris. This image has every dependency
# installed EXCEPT Ocellaris itself. Ocellaris can quickly be installed
# from the latest master branch branch by running
#
#   pip3 install --user git+https://bitbucket.org/trlandet/ocellaris.git
#
# See '.circleci/config.yml' for details on how Ocellaris is tested using
# this Docker image
#
# Reminder: how to build a new 'trlandet/fenics-dev:py3_CI' image:
# - Put this file in an empty directory. You can also use a non-empty
#   directory, but then it should not contain large files (like 
#   Singularity images etc.)
# - Build the Docker image
#       docker build .
# - Find the IMAGE_ID of the newly build image
#       docker images
# - Tag the image
#       docker tag IMAGE_ID trlandet/fenics-dev:py3_CI
# - Optionally test the image
#       docker run -it trlandet/fenics-dev:py3_CI
# - Login to Docker Hub
#       docker login -u trlandet
# - Push the updated image to Docker Hub
#       docker push trlandet/fenics-dev:py3_CI
FROM quay.io/fenicsproject/stable:2018.1.0.r3

USER root

# Install some prerequisites
RUN apt-get -qq update && \
    apt-get -y --with-new-pkgs \
    -o Dpkg::Options::="--force-confold" upgrade && \
    apt-get -y install python3-lxml && \
    apt-get clean && \
    update-alternatives --set mpirun /usr/bin/mpirun.mpich && \
    update-alternatives --set mpi /usr/include/mpich && \
    python3 -m pip install --upgrade pip && \
    python3 -m pip install --no-binary=h5py h5py && \
    python3 -m pip install PyYAML && \
    python3 -m pip install meshio && \
    python3 -m pip install yschema raschii && \
    cd && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
