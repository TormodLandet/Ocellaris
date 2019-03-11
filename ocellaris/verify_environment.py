# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import sys
import os
import warnings


def verify_env():
    """
    Check for dolfin and setup matplotlib backend
    """

    # Use non-GUI matplotlib backend if no GUI is available
    import matplotlib

    if has_tkinter():
        matplotlib.use('TkAgg')
    elif has_wx():
        matplotlib.use('WxAgg')
    else:
        matplotlib.use('Agg')

    error = False

    if not has_dolfin():
        msg = """

        ERROR: Could not import dolfin!
        Make sure FEniCS is properly installed
        Exiting due to error

        """
        sys.stderr.write(msg)
        error = True

    if not has_h5py():
        print('ERROR: missing h5py. Saving restart files will not work!', file=sys.stderr)
        error = True

    if not has_mpi4py():
        print('ERROR: missing mpi4py. FEniCS dolfin must be build with mpi4py', file=sys.stderr)
        error = True

    if not has_petsc4py():
        print('ERROR: missing petsc4py. FEniCS dolfin must be build with petsc4py', file=sys.stderr)
        error = True

    if not has_meshio():
        print('ERROR: missing meshio. Please install the "meshio" package!', file=sys.stderr)
        error = True

    if not has_yaml():
        print(
            'ERROR: missing required yaml module, please install the "PyYAML" package',
            file=sys.stderr,
        )
        error = True

    if error:
        exit()


def has_tkinter():
    try:
        import tkinter  # NOQA

        return True
    except ImportError:
        return False


def has_wx():
    try:
        import wx  # NOQA

        return True
    except Exception:
        return False


def has_yaml():
    try:
        import yaml  # NOQA

        return True
    except ImportError:
        return False


def has_dolfin():
    try:
        import dolfin  # NOQA

        return True
    except ImportError:
        return False


def has_meshio():
    try:
        import meshio  # NOQA

        return True
    except ImportError:
        return False


def has_mpi4py():
    try:
        import mpi4py  # NOQA

        return True
    except ImportError:
        return False


def has_petsc4py():
    try:
        from petsc4py import PETSc  # NOQA

        return True
    except ImportError:
        return False


def has_h5py():
    # Python cannot catch synchronous signals like SIGSEGV
    # so we need to warn before importing to avoid having a
    # mysterious crash at import time without any explaination

    # Warning
    msg1 = """
    ==========================================================
    WARNING: h5py is (probably) installed with an incompatible
    version of the hdf5 libraries. Make sure you install from
    source using the same HDF5 libs as dolfin.

    To reinstall h5py try something like this:

      pip3 install --no-binary=h5py h5py --user --upgrade \\
        --no-deps --ignore-installed

    Segfault will probably happen now!
    ==========================================================

    """

    # Apology
    msg2 = """
    ==========================================================
    Sorry about the h5py warning -- it seems that the HDF5
    libraries are compatible after all!

    Set the environment variable NOWARN_H5PY to stop this error
    message from being printed
    ==========================================================

    """

    # Check for h5py installation with bundled libraries
    probable_error = False
    for dirname in sys.path:
        hp = os.path.join(dirname, 'h5py')
        if os.path.isdir(os.path.join(hp, '.libs')):
            probable_error = True
        if os.path.isdir(hp):
            break

    # Do not warn if the user specifically wants us not to
    if os.getenv('NOWARN_H5PY'):
        probable_error = False

    # Warn that a SIGSEGV is forthcoming ...
    if probable_error:
        sys.stderr.write(msg1)
        sys.stderr.flush()

    try:
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=FutureWarning)
            import h5py  # NOQA
        return True
    except Exception:
        return False
    finally:
        if probable_error:
            sys.stderr.write(msg2)
