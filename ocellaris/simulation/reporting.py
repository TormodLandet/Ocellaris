from matplotlib import pyplot


class Reporting(object):
    def __init__(self, simulation):
        """
        Central place to register reports that will be output during
        the simulation
        """
        self.simulation = simulation
        self.timesteps = []
        self.timestep_xy_reports = {}
        simulation.hooks.add_pre_simulation_hook(self.setup_report_plotting, 'Reporting - setup plots')
    
    def setup_report_plotting(self):
        """
        Setup the reports to be shown in plots during the simulation
        """
        reps = self.simulation.input.get_value('reporting/reports_to_show', [], 'list(string)')
        self.figures = {}
        for report_name in reps:
            pyplot.ion()
            fig = pyplot.figure()
            ax = fig.add_subplot(111)
            ax.set_title(report_name)
            ax.set_xlabel('time')
            ax.set_ylabel(report_name)
            line, = ax.plot([], [])
            self.figures[report_name] = (fig, ax, line)
    
    def report_timestep_value(self, report_name, value):
        """
        Add a timestep to a report
        """
        time = self.simulation.time
        if not self.timesteps or not self.timesteps[-1] == time:
            self.timesteps.append(time)
        self.timestep_xy_reports.setdefault(report_name, []).append(value)
    
    def log_timestep_reports(self):
        """
        Write all reports for the finished time step to the log
        """
        info = []
        for report_name in sorted(self.timestep_xy_reports):
            value = self.timestep_xy_reports[report_name][-1]
            info.append('%s = %10g' % (report_name, value))
        it, t = self.simulation.timestep, self.simulation.time
        self.simulation.log.info('Reports for timestep = %5d, time = %10.4f, ' % (it, t) +
                                 ', '.join(info))
        
        for report_name in self.figures:
            fig, ax, line = self.figures[report_name]
            line.set_xdata(self.timesteps)
            line.set_ydata(self.timestep_xy_reports[report_name])
            ax.relim()
            ax.autoscale_view()
            fig.canvas.draw()
            fig.canvas.flush_events()
