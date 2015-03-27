
class Reporting(object):
    def __init__(self, simulation):
        """
        Central place to register reports that will be output during
        the simulation
        """
        self.simulation = simulation
        self.timesteps = []
        self.timestep_xy_reports = {}
    
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
