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
        print '\n    ERROR: Could not import dolfin!\n'
        print '    Make sure FEniCS is properly installed\n'
        print '    Exiting due to error\n'
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


def has_dolfin():
    try:
        import dolfin #@UnusedImport
        return True
    except ImportError:
        return False
