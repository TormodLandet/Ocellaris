import numpy
from matplotlib import pyplot
from . import Probe, register_probe

@register_probe('LineProbe')
class LineProbe(Probe):
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.input = probe_input

        # Read input
        name = self.input['name']
        startpos = self.input['startpos']
        endpos = self.input['endpos']
        N = self.input['Npoints']
        self.field_name = self.input['field']
        self.field = simulation.data[self.field_name]
        self.show_interval = self.input.get('show_interval', 0)
        self.show = self.show_interval != 0
        
        # Plot target values if provided
        self.has_target = False
        if 'target_abcissa' in self.input:
            self.has_target = True
            self.target_abcissa = self.input['target_abcissa']
            self.target_ordinate = self.input['target_ordinate']
            self.target_name = self.input.get('target_name', 'Target')
        
        # Handle 2D positions
        if len(startpos) == 2:
            startpos.append(0)
            endpos.append(0)
        
        # Get probe positions
        self.xvec = numpy.linspace(startpos[0], endpos[0], N)
        self.yvec = numpy.linspace(startpos[1], endpos[1], N)
        self.zvec = numpy.linspace(startpos[2], endpos[2], N)
        
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
        
        if not update_plot:
            return
        
        res = numpy.array([0.0])
        pos = numpy.array([0.0, 0.0, 0.0])
        
        probe_values = numpy.zeros_like(self.xvec)
        for i in range(len(probe_values)):
            pos[:] = (self.xvec[i], self.yvec[i], self.zvec[i])
            self.field.eval(res, pos)
            probe_values[i] = res[0]
        
        if self.xvec[0] != self.xvec[-1]:
            abcissa = self.xvec
            abcissa_label = 'x-pos'
        elif self.yvec[0] != self.yvec[-1]:
            abcissa = self.yvec
            abcissa_label = 'y-pos'
        else:
            abcissa = self.zvec
            abcissa_label = 'z-pos'
        
        if update_plot:
            self.line.set_xdata(abcissa)
            self.line.set_ydata(probe_values)
            self.ax.set_xlabel(abcissa_label)
            self.ax.relim()
            self.ax.autoscale_view()
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
