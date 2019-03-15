# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
from . import Probe, register_probe


@register_probe('LineProbe')
class LineProbe(Probe):
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.input = probe_input

        # Read input
        self.name = self.input.get_value('name', required_type='string')
        startpos = self.input.get_value('startpos', required_type='list(float)')
        endpos = self.input.get_value('endpos', required_type='list(float)')
        N = self.input.get_value('Npoints', required_type='int')
        self.field_name = self.input.get_value('field', required_type='string')
        self.field = simulation.data[self.field_name]

        # Should we write the data to a file
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        file_name = self.input.get_value('file_name', None, 'string')
        self.write_file = file_name is not None
        if self.write_file:
            self.file_name = prefix + '_' + file_name
            self.write_interval = self.input.get_value('write_interval', 1, 'int')

        # Handle 2D positions
        if len(startpos) == 2:
            startpos.append(0)
            endpos.append(0)

        # Get probe positions
        self.xvec = numpy.linspace(startpos[0], endpos[0], N)
        self.yvec = numpy.linspace(startpos[1], endpos[1], N)
        self.zvec = numpy.linspace(startpos[2], endpos[2], N)

        if self.write_file:
            self.output_file = open(self.file_name, 'wt')
            self.output_file.write('# Ocellaris line probe of the %s field\n' % self.field_name)
            self.output_file.write('# X = %s\n' % ' '.join('%15.5e' % x for x in self.xvec))
            self.output_file.write('# Y = %s\n' % ' '.join('%15.5e' % y for y in self.yvec))
            self.output_file.write('# Z = %s\n' % ' '.join('%15.5e' % z for z in self.zvec))
            self.output_file.write('#     time |-- probe values --> \n')

    def end_of_timestep(self):
        """
        Output the line probe at the end of the
        """
        it = self.simulation.timestep

        # Should we update the file?
        if not self.write_file or not (it == 1 or it % self.write_interval == 0):
            return

        # Get the value at the probe locations
        res = numpy.array([0.0])
        pos = numpy.array([0.0, 0.0, 0.0])
        probe_values = numpy.zeros_like(self.xvec)
        for i in range(len(probe_values)):
            pos[:] = (self.xvec[i], self.yvec[i], self.zvec[i])
            self.field.eval(res, pos)
            probe_values[i] = res[0]

        self.output_file.write(
            '%10.5f %s\n' % (self.simulation.time, ' '.join('%15.5e' % v for v in probe_values))
        )

    def end_of_simulation(self):
        """
        The simulation is done. Close the output file
        """
        if self.write_file:
            self.output_file.close()
