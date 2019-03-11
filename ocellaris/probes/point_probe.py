# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os
import dolfin
from ocellaris.utils import ocellaris_error, timeit
from . import Probe, register_probe


WRITE_INTERVAL = 1


@register_probe('PointProbe')
class PointProbe(Probe):
    description = 'Produce time series of point values'

    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.family = None
        self.degree = None

        # Read input
        inp = probe_input
        self.name = inp.get_value('name', required_type='string')
        N = len(inp.get_value('probe_points', required_type='list(any)'))
        self.probes = []
        for i in range(N):
            pconf = inp.get_value(
                'probe_points/%d' % i,
                required_type='list(any)',
                required_length=2 + simulation.ndim,
            )
            self.probes.append(pconf)

        # Verify input
        rank = simulation.rank
        probe_names = set()
        self.probes_on_process = []
        for pconf in self.probes:
            field_name = pconf[0]
            if not isinstance(field_name, str):
                ocellaris_error(
                    'Field name %r is not a string' % field_name,
                    'The first point probe parameter must be a string',
                )
            if field_name not in simulation.data:
                ocellaris_error(
                    'Field %r does not exist' % field_name,
                    'The first point probe parameter must be the name of a field',
                )
            pname = pconf[1]
            if not isinstance(pname, str):
                ocellaris_error(
                    'Probe name %r is not a string' % pname,
                    'The second point probe parameter must be a string',
                )
            if pname in probe_names:
                ocellaris_error(
                    'Probe name %r is duplicate' % pname,
                    'The second point probe parameter must be a unique string',
                )

            # Check if this point is on this process
            coords = tuple(pconf[2:])
            simulation.log.info('        Probe %s of %s at %r' % (pname, field_name, coords))
            pt = dolfin.Point(*coords)
            field = simulation.data[field_name]
            try:
                field(pt)
                self.probes_on_process.append((field_name, pname, pt))
                simulation.log.info('            Probe %s found on rank %r' % (pname, rank))
            except Exception as e:
                simulation.log.info('            Probe %s not found on rank %r' % (pname, rank))

        # File name to write output to
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        self.file_name = inp.get_value('file_name', prefix + '_%s.tsv' % self.name, 'string')

        if rank == 0:
            if os.path.isfile(self.file_name):
                simulation.log.info('        Appending to TSV file %s' % self.file_name)
                self.tsv = open(self.file_name, 'a')
            else:
                simulation.log.info('        Creating TSV file %s' % self.file_name)
                self.tsv = open(self.file_name, 'w')
                head = ['t'] + [pconf[1] for pconf in self.probes]
                self.tsv.write('\t'.join(head) + '\n')
        else:
            self.tsv = None

        # Add field to list of IO plotters and listen for flush events
        inp_key = probe_input.basepath + 'write_interval'
        simulation.io.add_plotter(self.write_field, inp_key, WRITE_INTERVAL)
        self.simulation.hooks.add_custom_hook('flush', self.flush, 'Flush tsv file')

    @timeit.named('PointProbe.write_field')
    def write_field(self):
        """
        Find and output the plane probe
        """
        sim = self.simulation

        # Query values on this process
        rank_values = []
        for field_name, probe_name, point in self.probes_on_process:
            field = sim.data[field_name]
            field_val = field(point)
            rank_values.append((probe_name, field_val))

        # Send all values to the root process
        comm = sim.data['mesh'].mpi_comm()
        all_values = comm.gather(rank_values)

        # Append to output file
        if all_values is not None:
            # Get all distributed values
            values = {}
            for rank_vals in all_values:
                for k, v in rank_vals:
                    if k in values and v != values[k]:
                        sim.log.warning(
                            'Probe %s got different values %r and %r' % (k, v, values[k])
                        )
                    values[k] = v

            # Construct data list and write to file
            data = [repr(sim.time)]
            for pconf in self.probes:
                probe_name = pconf[1]
                if probe_name not in values:
                    sim.log.warning(
                        'Probe %s not found on any process! Writing dummy value -1e100!'
                        % probe_name
                    )
                    data.append(repr(-1e100))
                else:
                    data.append(repr(values[probe_name]))
            self.tsv.write('\t'.join(data) + '\n')

    def flush(self):
        if self.tsv is not None:
            self.tsv.flush()
