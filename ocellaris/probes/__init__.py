# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import time
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
        ocellaris_error(
            'Postprocessing probe "%s" not found' % name,
            'Available probe:\n'
            + '\n'.join(
                '  %-20s - %s' % (n, s.description) for n, s in sorted(_PROBES.items())
            ),
        )
        raise


def setup_probes(simulation):
    """
    Install probes from a simulation input
    """

    def hook_start(probe):
        return lambda: probe.new_simulation()

    def hook_pre_timestep(probe):
        return lambda timestep_number, t, dt: probe.new_timestep(timestep_number, t, dt)

    def hook_post_timestep(probe):
        return lambda: probe.end_of_timestep()

    def hook_final(probe):
        return lambda success: probe.end_of_simulation()

    simulation.probes = {}
    Nprobes = len(simulation.input.get_value('probes', []))

    if Nprobes:
        simulation.log.info('Setting up probes')

    for i in range(Nprobes):
        # Read input about this probe
        inp = simulation.input.get_value('probes/%d' % i, required_type='Input')
        probe_name = inp.get_value('name', 'unnamed', 'string')
        probe_type = inp.get_value('type', required_type='string')
        enabled = inp.get_value('enabled', True, 'bool')

        if not enabled:
            continue
        t1 = time.time()
        simulation.log.info('    Setting up %s probe %s' % (probe_type, probe_name))

        # Get the probe object
        probe_class = get_probe(probe_type)
        probe = probe_class(simulation, inp)

        # Store for access from other parts of the code
        if probe_name in simulation.probes:
            simulation.log.warning(
                '    Probe name "%s" used for multiple probes. Names should be unique!'
            )
        simulation.probes[probe_name] = probe

        # Register callbacks
        if probe.new_simulation is not None:
            simulation.hooks.add_pre_simulation_hook(
                hook_start(probe), 'Probe "%s"' % probe_name
            )
        if probe.new_timestep is not None:
            simulation.hooks.add_pre_timestep_hook(
                hook_pre_timestep(probe), 'Probe "%s"' % probe_name
            )
        if probe.end_of_timestep is not None:
            simulation.hooks.add_post_timestep_hook(
                hook_post_timestep(probe), 'Probe "%s"' % probe_name
            )
        if probe.end_of_simulation is not None:
            simulation.hooks.add_post_simulation_hook(
                hook_final(probe), 'Probe "%s"' % probe_name
            )

        simulation.log.info(
            '    Probe %s set-up done in %.1f seconds' % (probe_name, time.time() - t1)
        )


class Probe(object):
    description = 'No description available'

    def __init__(self, simulation, probe_input):
        """
        A base class for post-processing probes
        """
        self.simulation = simulation
        self.input = probe_input

    # These can be callables (methods) in derived classes
    new_simulation = None
    new_timestep = None
    end_of_timestep = None
    end_of_simulation = None


from . import line_probe
from . import iso_surface
from . import plane_probe
from . import point_probe
