from ocellaris.utils import ocellaris_error

_PROBES = {}

def add_probe(name, probe_class):
    """
    Register a postprocessing probe
    """
    _PROBES[name] = probe_class

def register_probe(name):
    """
    A class decorator to register postprocessing probes
    """
    def register(probe_class):
        add_probe(name, probe_class)
        return probe_class
    return register

def get_probe(name):
    """
    Return a postprocessing probe by name
    """
    try:
        return _PROBES[name]
    except KeyError:
        ocellaris_error('Postprocessing probe "%s" not found' % name,
                        'Available probe:\n' +
                        '\n'.join('  %-20s - %s' % (n, s.description) 
                               for n, s in sorted(_PROBES.items())))
        raise

def setup_probes(simulation):
    """
    Install probes from a simulation input
    """
    def hook_timestep(probe):
        return lambda: probe.end_of_timestep()
    
    def hook_final(probe):
        return lambda success: probe.end_of_simulation()
    
    Nprobes = len(simulation.input.get_value('probes', []))
    for i in range(Nprobes):
        inp = simulation.input.get_value('probes/%d' % i, required_type='Input')
        probe_name = inp.get_value('name', 'unnamed', 'string')
        probe_type = inp.get_value('type', required_type='string')
        probe_class = get_probe(probe_type)
        probe = probe_class(simulation, inp)
        simulation.hooks.add_post_timestep_hook(hook_timestep(probe), 'Probe "%s"' % probe_name)
        simulation.hooks.add_post_simulation_hook(hook_final(probe), 'Probe "%s"' % probe_name)

class Probe(object):
    def __init__(self, simulation, probe_input):
        """
        A base class for post-processing probes
        """
        self.simulation = simulation
        self.input = probe_input
    
    def end_of_timestep(self):
        pass
    
    def end_of_simulation(self, success):
        pass

from . import line_probe
from . import iso_surface
