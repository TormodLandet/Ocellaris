# Check for presence of FEniCS etc
from .verify_environment import verify_env
verify_env()


__version__ = '2018.1.0.dev0'


# This should potentially be made local to the mesh creation routines
import dolfin
dolfin.parameters['ghost_mode'] = 'shared_vertex'
del dolfin


def get_version():
    """
    Return the version number of Ocellaris
    """
    return __version__


def get_detailed_version():
    """
    Return the version number of Ocellaris including
    source control commit revision information
    """
    import os, subprocess
    this_dir = os.path.dirname(os.path.abspath(__file__))
    proj_dir = os.path.abspath(os.path.join(this_dir, '..'))
    if os.path.isdir(os.path.join(proj_dir, '.git')):
        cmd = ['git', 'describe', '--always']
        version = subprocess.check_output(cmd, cwd=proj_dir)
        local_version = '+git.' + version.decode('utf8').strip()
    else:
        local_version = ''
    return get_version() + local_version


# Convenience imports for scripting
from .simulation import Simulation
from .run import setup_simulation, run_simulation
