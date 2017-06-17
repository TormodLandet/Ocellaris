Bootstrap: docker
From: trlandet/fenics-dev:latest

%post
    # Install additional dependencies (FEniCS is allready installed)
    pip2 install PyYAML
    pip3 install PyYAML

    # Install the solenoidal limiter
    # (this should get pip installable if it becomes a permanent dependency ...)
    cd /opt/src
    git clone https://bitbucket.org/trlandet/divergencefreelimiter.git
    cp -R divergencefreelimiter/solenoidal /usr/local/lib/python2.7/dist-packages

    # Install Ocellaris
    cd /opt/src
    git clone https://bitbucket.org/trlandet/ocellaris.git
    cd ocellaris
    pip2 install -e .

    # Set matplotlib default backend to something that works
    sed -i -e 's/backend[[:space:]]*:[[:space:]]TkAgg/backend: Agg/g' /usr/local/lib/python2.7/dist-packages/matplotlib/mpl-data/matplotlibrc
    sed -i -e 's/backend[[:space:]]*:[[:space:]]TkAgg/backend: Agg/g' /usr/local/lib/python3.5/dist-packages/matplotlib/mpl-data/matplotlibrc

    # General setup
    chmod a+rwX -R /opt
    echo 'PS1="Ocellaris \w> "' > /opt/bashrc

%runscript
    echo "Welcome to the Ocellaris Singularity container"
    echo "FEniCS is installed for python2 and python3"
    echo "Ocellaris is installed for python2"
    export FENICS_PREFIX=${FENICS_PREFIX:-/opt}
    export FENICS_SRC_DIR=${FENICS_PREFIX}/src
    export SLEPC_DIR=${FENICS_PREFIX}
    export PETSC_DIR=${FENICS_PREFIX}
    export PATH=${FENICS_PREFIX}/bin:/usr/lib/ccache:$PATH
    export LD_LIBRARY_PATH=${FENICS_PREFIX}/lib:$LD_LIBRARY_PATH
    export MANPATH=${FENICS_PREFIX}/share/man:$MANPATH
    export CPATH=${FENICS_PREFIX}/include:$CPATH
    touch ~/.sudo_as_admin_successful
    exec /bin/bash --rcfile /opt/bashrc "$@"

