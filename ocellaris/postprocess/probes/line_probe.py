import numpy
from matplotlib import pyplot
from . import Probe, register_probe


@register_probe('LineProbe')
class LineProbe(Probe):
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.input = probe_input

        # Read input
        name = self.input.get_value('name', required_type='string')
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
        
        # Should we pop up a matplotlib window when running?
        self.show_interval = self.input.get_value('show_interval', 0, 'float')
        self.show = self.show_interval != 0
        
        # Plot target values if provided
        self.has_target = False
        if 'target_abcissa' in self.input:
            self.has_target = True
            self.target_abcissa = self.input.get_value('target_abcissa', required_type='list(float)')
            self.target_ordinate = self.input.get_value('target_ordinate', required_type='list(float)')
            self.target_name = self.input.get_value('target_name', 'Target', 'string')
        
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
        
        if self.show:
            pyplot.ion()
            self.fig = pyplot.figure()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_title('Line probe %s' % name)
            self.ax.set_ylabel(self.field_name)
            self.line, = self.ax.plot([], [])
            
            if self.has_target:
                self.target_line, = self.ax.plot(self.target_abcissa, self.target_ordinate, 'kv')
                self.ax.legend(['Ocellaris', self.target_name], loc='best')
        
    def end_of_timestep(self):
        """
        Output the line probe at the end of the
        """
        it = self.simulation.timestep
        
        # Should we update the plot?
        update_plot = False 
        if self.show and (it == 1 or it % self.show_interval == 0):
            update_plot = True

        # Should we update the file?
        update_file = False
        if self.write_file and (it == 1 or it % self.write_interval == 0):
            update_file = True
        
        # Do not do any postprocessing for non-requested time steps
        if not (update_file or update_plot):
            return
        
        # Get the value at the probe locations
        res = numpy.array([0.0])
        pos = numpy.array([0.0, 0.0, 0.0])
        probe_values = numpy.zeros_like(self.xvec)
        for i in range(len(probe_values)):
            pos[:] = (self.xvec[i], self.yvec[i], self.zvec[i])
            self.field.eval(res, pos)
            probe_values[i] = res[0]
        
        # For plotting, figure out which axis is the abcissa
        if self.xvec[0] != self.xvec[-1]:
            abcissa = self.xvec
            abcissa_label = 'x-pos'
        elif self.yvec[0] != self.yvec[-1]:
            abcissa = self.yvec
            abcissa_label = 'y-pos'
        else:
            abcissa = self.zvec
            abcissa_label = 'z-pos'
        
        if update_file:
            self.output_file.write('%10.5f %s\n' % (self.simulation.time,
                                   ' '.join('%15.5e' % v for v in probe_values)))
        
        if update_plot:
            self.line.set_data(abcissa, probe_values)
            self.ax.set_xlabel(abcissa_label)
            self.ax.relim()
            self.ax.autoscale_view()
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
    
    def end_of_simulation(self):
        """
        The simulation is done. Close the output file
        """
        if self.write_file:
            self.output_file.close()
