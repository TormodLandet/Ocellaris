import json
import collections
import time
import dolfin
from .plot import Plotter
from .utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data

class Simulation(object):
    def __init__(self):
        self.input = {}
        self.data = {}
        self.plotting = Plotting(self)
        self.reporting = Reporting(self)
        self.log = Log(self)
        self._pre_timestep_hooks = []
        self._post_timestep_hooks = []
        
        # Several parts of the code wants to know these things,
        # so we keep them in a central place
        self.ndim = 0
        self.timestep = 0
        self.time = 0.0
        self.dt = 0.0
        
        # For timing the analysis
        self.prevtime = self.starttime = time.time()
        
        # Setup reporting
        rep = lambda report: self.reporting.print_timestep_report() if report else None
        self.add_post_timestep_hook(rep)
    
    def set_mesh(self, mesh):
        """
        Set the computational domain
        """
        self.data['mesh'] = mesh
        self.ndim = mesh.topology().dim()
        self.update_mesh_data()
    
    def update_mesh_data(self):
        """
        Some precomputed values must be calculated before the timestepping
        and updated every time the mesh changes
        """
        init_connectivity(self)
        precompute_cell_data(self)
        precompute_facet_data(self)
    
    def add_pre_timestep_hook(self, hook):
        """
        Add a function that will run before the solver in each time step
        """
        self._pre_timestep_hooks.append(hook)
        
    def add_post_timestep_hook(self, hook):
        """
        Add a function that will run after the solver in each time step
        """
        self._post_timestep_hooks.append(hook)
    
    def new_timestep(self, timestep_number, t, dt):
        """
        Called at the start of a new time step
        """
        self.timestep = timestep_number
        self.time = t
        self.dt = dt
        for hook in self._pre_timestep_hooks:
            hook(t, dt)
    
    def end_timestep(self, report=True):
        """
        Called at the end of a time step
        """
        # Report the time spent in this time step
        newtime = time.time()
        self.reporting.report_timestep_value('tstime', newtime-self.prevtime)
        self.reporting.report_timestep_value('tottime', newtime-self.starttime)
        self.prevtime = newtime
        
        for hook in self._post_timestep_hooks:
            hook(report)
    
    def read_json_input_file(self, filename):
        with open(filename, 'rt') as inpf:
            inp = json.load(inpf, object_pairs_hook=collections.OrderedDict)
        assert inp['program'] == 'ocellaris'
        assert inp['version'] == 1.0
        assert inp['type'] == 'input'
        self.input = inp

class Plotting(object):
    def __init__(self, simulation):
        """
        Central place to register plots that will be output during
        the simulation
        """
        self.simulation = simulation
        self.plots = {}
    
    def add_plot(self, plot_name, plotter):
        """
        Add a plot to the simulation
        """
        if not hasattr(plotter, 'plot'):
            # This is not a plotter but something that can be plotted
            plotter = Plotter(plotter)
        
        self.plots[plot_name] = plotter
    
    def plot_all(self):
        """
        Plot all registered plotters
        """
        for name in self.plots:
            self.plot(name)
    
    def plot(self, name, extra=''):
        """
        Plot the plotter with the given name
        """
        sim = self.simulation
        name_template_dafault = 'fig/{name}_{timestep:07d}_{t:010.6f}{extra}.png'
        name_template = sim.input.get('plots', {}).get('name_template', name_template_dafault)
        filename = name_template.format(name=name, timestep=sim.timestep, t=sim.time, extra=extra)
        self.plots[name].plot(filename)
        
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
    
    def print_timestep_report(self):
        """
        Plot all registered plotters
        """
        info = []
        for report_name in sorted(self.timestep_xy_reports):
            value = self.timestep_xy_reports[report_name][-1]
            info.append('%s = %10g' % (report_name, value))
        it, t = self.simulation.timestep, self.simulation.time
        print 'Reports for timestep = %5d, time = %10.4f,' % (it, t),
        print ', '.join(info)

class Log(object):
    def __init__(self, simulation):
        self.simulation = simulation
        self.log_level = dolfin.INFO
    
    def set_log_level(self, log_level):
        """
        Set the Ocillaris log level
        (not the dolfin log level!)
        """
        self.log_level = log_level
    
    def info(self, text):
        "Log info"
        if self.log_level <= dolfin.INFO:
            print text
    
    def debug(self, text):
        "Log debug"
        if self.log_level <= dolfin.DEBUG:
            print text
