import sys, os

def verify_env():
    """
    Check for dolfin and setup matplotlib backend
    """
    
    # Use non-GUI matplotlib backend if no GUI is available
    import matplotlib
    if has_tk():
        matplotlib.use('TkAgg')
    elif has_wx():
        matplotlib.use('WxAgg')
    else:
        matplotlib.use('Agg')
    
    if not has_dolfin():
        msg = """
        
        ERROR: Could not import dolfin!
        Make sure FEniCS is properly installed
        Exiting due to error
        
        """
        sys.stderr.write(msg)
        exit()
    
    if not has_h5py():
        print('WARNING: missing h5py. Saving restart files will not work!')
    
    if not has_yaml():
        print('Missing required yaml module, please install the PyYAML package',
              file=sys.stderr)    
        exit()


def has_tk():
    try:
        from six.moves import tkinter #@UnusedImport
        return True
    except ImportError:
        return False


def has_wx():
    try:
        import wx #@UnusedImport
        return True
    except ImportError:
        return False


def has_yaml():
    try:
        import yaml #@UnusedImport
        return True
    except ImportError:
        return False

def has_dolfin():
    try:
        import dolfin #@UnusedImport
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
        import h5py #@UnusedImport
        return True
    except:
        return False
    finally:
        if probable_error:
            sys.stderr.write(msg2)
