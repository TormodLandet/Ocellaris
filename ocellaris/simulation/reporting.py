# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from collections import OrderedDict
from matplotlib import pyplot
from ocellaris.utils import ocellaris_error, timeit


class Reporting(object):
    def __init__(self, simulation):
        """
        Central place to register reports that will be output during
        the simulation
        """
        self.simulation = simulation
        self.timesteps = []
        self.timestep_xy_reports = OrderedDict()
        simulation.hooks.add_pre_simulation_hook(
            self.setup_report_plotting, 'Reporting - setup plots'
        )

    def setup_report_plotting(self):
        """
        Setup the reports to be shown in plots during the simulation
        """
        if self.simulation.rank != 0:
            return  # Do not plot on non root processes

        reps = self.simulation.input.get_value(
            'reporting/reports_to_show', [], 'list(string)'
        )
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
        rep = self.timestep_xy_reports.setdefault(report_name, [])
        rep.extend([None] * (len(self.timesteps) - len(rep)))
        rep[-1] = value

    def get_report(self, report_name):
        """
        Get a the time series of a reported value
        """
        t = self.timesteps
        rep = self.timestep_xy_reports[report_name]
        return t[: len(rep)], rep

    @timeit.named('reporting log_timestep_reports')
    def log_timestep_reports(self):
        """
        Write all reports for the finished time step to the log
        """
        info = []
        for report_name in self.timestep_xy_reports:
            value = self.timestep_xy_reports[report_name][-1]
            info.append('%s = %10g' % (report_name, value))
        it, t = self.simulation.timestep, self.simulation.time
        self.simulation.log.info(
            'Reports for timestep = %5d, time = %10.4f, ' % (it, t) + ', '.join(info)
        )

        # Update interactive report plots
        self._update_plots()

    def _update_plots(self):
        """
        Update plots requested in input (reporting/reports_to_show)
        """
        if self.simulation.rank != 0:
            return  # Do not plot on non root processes

        for report_name in self.figures:
            if not report_name in self.timestep_xy_reports:
                ocellaris_error(
                    'Unknown report name: "%s"' % report_name,
                    'Cannot plot this report, it does not exist',
                )

            abscissa = self.timesteps
            ordinate = self.timestep_xy_reports[report_name]

            if len(abscissa) != len(ordinate):
                ocellaris_error(
                    'Malformed report data: "%s"' % report_name,
                    'Length of t is %d while report is %d'
                    % (len(abscissa), len(ordinate)),
                )

            fig, ax, line = self.figures[report_name]
            line.set_xdata(abscissa)
            line.set_ydata(ordinate)
            ax.relim()
            ax.autoscale_view()
            fig.canvas.draw()
            fig.canvas.flush_events()
