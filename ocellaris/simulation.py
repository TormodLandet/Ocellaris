import json
import collections
import time
import dolfin
from .postprocess import Plotter
from .utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data
from .utils import timeit, report_error
from sympy.functions.combinatorial import numbers

class Simulation(object):
    def __init__(self):
        """
        Represents one Navier-Stokes simulation. The simulation 
        connects the input file, geometry, mesh with the solver
        and the result plotting and reporting tools     
        """
        self.input = Input(self)
        self.data = {}        
        self.plotting = Plotting(self)
        self.reporting = Reporting(self)
        self.log = Log(self)
        self._pre_timestep_hooks = []
        self._post_timestep_hooks = []
        self._post_simulation_hooks = []
        
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
    
    def add_post_simulation_hook(self, hook):
        """
        Add a function that will run after the simulation is done
        """
        self._post_simulation_hooks.append(hook)
    
    @timeit
    def new_timestep(self, timestep_number, t, dt):
        """
        Called at the start of a new time step
        """
        self.timestep = timestep_number
        self.time = t
        self.dt = dt
        for hook in self._pre_timestep_hooks:
            hook(t, dt)
    
    @timeit
    def end_timestep(self, report=True):
        """
        Called at the end of a time step
        
        Arguments:
            report: True if something should be written to
                console to summarise the last time step
        """
        # Report the time spent in this time step
        newtime = time.time()
        self.reporting.report_timestep_value('tstime', newtime-self.prevtime)
        self.reporting.report_timestep_value('tottime', newtime-self.starttime)
        self.prevtime = newtime
        
        # Report the maximum velocity
        vels = 0
        for d in range(self.ndim):
            vels += self.data['u'][d].vector().array()**2
        vel_max = vels.max()**0.5
        self.reporting.report_timestep_value('umax', vel_max)
        
        for hook in self._post_timestep_hooks:
            hook(report)
    
    def end_simulation(self, success):
        """
        Report that  the simulation is done.
        
        Arguments:
            success: True if nothing went wrong, False for
            diverging solution and other problems
        """
        for hook in self._post_simulation_hooks:
            hook(success)

    
class Input(dict):
    UNDEFINED = object()
    
    def __init__(self, simulation):
        """
        Holds the input values provided by the user
        """
        self.simulation = simulation
    
    def read_json(self, filename):
        """
        Read the input to an Ocellaris simulation from a JSON 
        formated input file. The user will get an error if the
        input file is malformed 
        """
        try:
            with open(filename, 'rt') as inpf:
                inp = json.load(inpf, object_pairs_hook=collections.OrderedDict)
        except ValueError as e:
            report_error('Error on input file', str(e))
        
        assert inp['program'] == 'ocellaris'
        assert inp['version'] == 1.0
        assert inp['type'] == 'input'
        self.clear()
        self.update(inp)
        
    def get_value(self, path, default_value=UNDEFINED, required_type='any'):
        """
        Get an input value by its path in the input dictionary
        
        Gives an error if there is no default value supplied
        and the  input variable does not exist
        
        Arguments:
            path: a list of path components or the "/" separated
                path to the variable in the input dictionary
            default_value: the value to return if the path does
                not exist in the input dictionary
            required_type: expected type of the variable. Giving 
                type="any" does no type checking
        
        Returns:
            The input value if it exist otherwise the default value
        """
        # Allow path to be a list or a "/" separated string
        if isinstance(path, basestring):
            pathstr = path
            path = path.split('/')
        else:
            pathstr = '/'.join(path)
        
        d = self
        for p in path:
            if p not in d:
                if default_value is self.UNDEFINED:
                    report_error('Missing parameter on input file',
                                 'Missing required input parameter:\n  %s' % pathstr)
                else:
                    return default_value
            d = d[p]
        
        def check_isinstance(value, classes):
            """
            Give error if the input data is not of the required type
            """
            if not isinstance(value, classes):
                report_error('Malformed data on input file',
                             'Parameter %s should be of type %s,\nfound %r %r' % 
                             (pathstr, required_type, value, type(value)))
        
        # Validate according to required data type
        number = (int, long, float)
        if required_type == 'bool':
            check_isinstance(d, bool)
        elif required_type == 'float':
            check_isinstance(d, number)
        elif required_type == 'int':
            check_isinstance(d, int)
        elif required_type == 'string':
            check_isinstance(d, basestring)
        elif required_type == 'list(float)':
            check_isinstance(d, list)
            for elem in d:
                check_isinstance(elem, number)
        elif required_type == 'list(dict)':
            check_isinstance(d, list)
            for elem in d:
                check_isinstance(elem, dict)
        elif required_type == 'any':
            pass
        else:
            raise ValueError('Unknown required_type %s' % required_type)
        return d
    
    def get_output_path(self, path, default_value=UNDEFINED):
        """
        Get the name of an output file
        """
        prefix = self.get_value('output/prefix', '')
        filename = self.get_input(path, default_value)
        return prefix + filename

class Plotting(object):
    def __init__(self, simulation):
        """
        Central place to register plots that will be output during
        the simulation
        """
        self.simulation = simulation
        self.plots = {}
    
    def add_plot(self, plot_name, plotter, **options):
        """
        Add a plot to the simulation
        """
        if not hasattr(plotter, 'plot'):
            # This is not a plotter but something that can be plotted
            plotter = Plotter(self.simulation, plotter, **options)
        
        self.plots[plot_name] = plotter
        return plotter
    
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
        name_template = sim.input.get_value('plots/name_template', name_template_dafault, 'string')
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
        Set the Ocellaris log level
        (not the dolfin log level!)
        """
        self.log_level = log_level
        
    def warning(self, text):
        "Log warning"
        if self.log_level <= dolfin.WARNING:
            print text
    
    def info(self, text):
        "Log info"
        if self.log_level <= dolfin.INFO:
            print text
    
    def debug(self, text):
        "Log debug"
        if self.log_level <= dolfin.DEBUG:
            print text
