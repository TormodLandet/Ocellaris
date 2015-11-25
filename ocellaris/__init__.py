try:
    import dolfin #@UnusedImport
except ImportError:
    print '\n    ERROR: Could not import dolfin!\n'
    print '    Make sure FEniCS is properly installed\n'
    print '    Exiting due to error\n'
    exit()

_VERSION = '0.1'
def get_version():
    """
    Return the version number of Ocellaris
    """
    return _VERSION

def get_detailed_version():
    """
    Return the version number of Ocellaris including
    source control commit revision information
    """
    import os, subprocess
    this_dir = os.path.dirname(os.path.abspath(__file__))
    proj_dir = os.path.abspath(os.path.join(this_dir, '..'))
    if os.path.isdir(os.path.join(proj_dir, '.hg')):
        cmd = ['hg', 'log', '-r', '.', '--template', 
               '{latesttag}-{latesttagdistance}-{node|short}']
        version = subprocess.check_output(cmd)
        return version.strip()

# Convenience imports for scripting
from .postprocess import Plotter
from .simulation import Simulation
from .run import setup_simulation, run_simulation
