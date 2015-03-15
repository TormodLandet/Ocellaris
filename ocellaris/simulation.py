import os, collections, time
import yaml
import dolfin
from .postprocess import Plotter
from .utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data
from .utils import timeit, report_error

class Simulation(object):
    def __init__(self):
        """
        Represents one Ocellaris simulation. The Simulation class 
        connects the input file, geometry, mesh and more with the
        solver and the result plotting and reporting tools     
        """
        self.hooks = Hooks(self)
        self.input = Input(self)
        self.data = {}        
        self.plotting = Plotting(self)
        self.reporting = Reporting(self)
        self.log = Log(self)
        
        # Several parts of the code wants to know these things,
        # so we keep them in a central place
        self.ndim = 0
        self.timestep = 0
        self.time = 0.0
        self.dt = 0.0
        
        # These will be filled out when ocellaris.run is setting up
        # the solver. Included here for documentation purposes only
        self.solver = None
        self.multi_phase_model = None
        
        # For timing the analysis
        self.prevtime = self.starttime = time.time()
        self.hooks.add_pre_timestep_hook(self._at_start_of_timestep)
        self.hooks.add_post_timestep_hook(self._at_end_of_timestep)
    
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
        
    def _at_start_of_timestep(self, timestep_number, t, dt):
        self.timestep = timestep_number
        self.time = t
        self.dt = dt
    
    def _at_end_of_timestep(self, report):
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

class Hooks(object):
    def __init__(self, simulation):
        """
        This class allows registering functions to run at
        given times during the simulation, e.g. to update
        some values for the next time step, report something
        after each time step or clean up after the simulation
        """
        self.simulation = simulation
        self._pre_simulation_hooks = []
        self._pre_timestep_hooks = []
        self._post_timestep_hooks = []
        self._post_simulation_hooks = []
    
    def add_pre_simulation_hook(self, hook):
        """
        Add a function that will run before the simulation starts
        """
        self._pre_simulation_hooks.append(hook)
    
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
        
    def simulation_started(self):
        """
        Called by the solver when the simulation starts
        
        Will run all pre simulation hooks in the reverse
        order they have been added
        """
        for hook in self._pre_simulation_hooks[::-1]:
            hook()
    
    @timeit
    def new_timestep(self, timestep_number, t, dt):
        """
        Called by the solver at the beginning of a new time step
        
        Will run all pre timestep hooks in the reverse
        order they have been added 
        """
        for hook in self._pre_timestep_hooks[::-1]:
            hook(timestep_number, t, dt)
    
    @timeit
    def end_timestep(self):
        """
        Called by the solver at the end of a time step
        
        Will run all post timestep hooks in the reverse
        order they have been added
        """
        report = True
        for hook in self._post_timestep_hooks[::-1]:
            hook(report)
    
    def simulation_ended(self, success):
        """
        Called by the solver when the simulation is done
        
        Will run all post simulation hooks in the reverse
        order they have been added
        
        Arguments:
            success: True if nothing went wrong, False for
            diverging solution and other problems
        """
        self.simulation.success = success
        for hook in self._post_simulation_hooks[::-1]:
            hook(success)

class UndefinedParameter(object):
    def __repr__(self):
        "For Sphinx"
        return '<UNDEFINED>'
UNDEFINED = UndefinedParameter()


class Input(collections.OrderedDict):
    def __init__(self, simulation):
        """
        Holds the input values provided by the user
        """
        super(Input, self).__init__()
        self.simulation = simulation
    
    def read_yaml(self, file_name):
        """
        Read the input to an Ocellaris simulation from a YAML 
        formated input file. The user will get an error if the
        input file is malformed 
        """
        self._setup_yaml()
        try:
            with open(file_name, 'rt') as inpf:
                inp = yaml.load(inpf)
        except ValueError as e:
            report_error('Error on input file', str(e))
        
        assert 'ocellaris' in inp
        assert inp['ocellaris']['type'] == 'input'
        assert inp['ocellaris']['version'] == 1.0
        
        self.clear()
        self.update(inp)
        self.file_name = file_name
    
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
                if default_value is UNDEFINED:
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
        dict_types = (dict, collections.OrderedDict)
        if required_type == 'bool':
            check_isinstance(d, bool)
        elif required_type == 'float':
            check_isinstance(d, number)
        elif required_type == 'int':
            check_isinstance(d, int)
        elif required_type == 'string':
            check_isinstance(d, basestring)
            # SWIG does not like Python 2 Unicode objects
            d = str(d)
        elif required_type == 'dict(string:dict)':
            check_isinstance(d, dict_types)
            for key, val in d.items():
                check_isinstance(key, basestring)
                check_isinstance(val, dict_types)
        elif required_type == 'dict(string:list)':
            check_isinstance(d, dict_types)
            for key, val in d.items():
                check_isinstance(key, basestring)
                check_isinstance(val, list)
        elif required_type == 'list(float)':
            check_isinstance(d, list)
            for elem in d:
                check_isinstance(elem, number)
        elif required_type == 'list(dict)':
            check_isinstance(d, list)
            for elem in d:
                check_isinstance(elem, dict_types)
        elif required_type == 'any':
            pass
        else:
            raise ValueError('Unknown required_type %s' % required_type)
        return d
    
    def get_output_file_path(self, path, default_value=UNDEFINED):
        """
        Get the name of an output file
        
        Automatically prefixes the file name with the output prefix
        """
        prefix = self.get_value('output/prefix', '')
        filename = self.get_value(path, default_value, 'string')
        if default_value is None and filename is None:
            return None
        else:
            return prefix + filename
        
    def get_input_file_path(self, file_name):
        """
        Serch first relative to the current working dir and then
        relative to the input file dir
        """
        # Check if the path is absolute or relative to the
        # working directory
        if os.path.exists(file_name):
            return file_name
        self.simulation.log.debug('File does not exist: %s' % file_name)
        
        # Check if the path is relative to the inouf file dir
        inp_file_dir = os.path.dirname(self.file_name)
        pth2 = os.path.join(inp_file_dir, file_name)
        if os.path.exists(pth2):
            return pth2
        self.simulation.log.debug('File does not exist: %s' % pth2)
        
        report_error('File not found', 'The specified file "%s" was not found' % file_name)
    
    def _setup_yaml(self):
        """
        Make PyYaml load and store keys in dictionaries 
        ordered like they were on the input file
        """
        _mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
    
        def dict_representer(dumper, data):
            return dumper.represent_dict(data.iteritems())
    
        def dict_constructor(loader, node):
            return collections.OrderedDict(loader.construct_pairs(node))
    
        yaml.add_representer(collections.OrderedDict, dict_representer)
        yaml.add_constructor(_mapping_tag, dict_constructor)


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
        
        # Setup reporting after each time step
        rep = lambda report: self.log_timestep_reports() if report else None
        self.simulation.hooks.add_post_timestep_hook(rep)
    
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


class Log(object):
    def __init__(self, simulation):
        self.simulation = simulation
        self.log_level = dolfin.INFO
        self.simulation.hooks.add_post_simulation_hook(lambda success: self.end_of_simulation())
        self.write_log = False
    
    def _write(self, message):
        """
        Internal helper to write message to
        console and log file
        """
        if self.write_log:
            self.log_file.write(message + '\n')
        print message
    
    def set_log_level(self, log_level):
        """
        Set the Ocellaris log level
        (not the dolfin log level!)
        """
        self.log_level = log_level
    
    def error(self, message):
        "Log an error message"
        if self.log_level <= dolfin.ERROR:
            self._write(message)
    
    def warning(self, message=''):
        "Log a warning message"
        if self.log_level <= dolfin.WARNING:
            self._write(message)
    
    def info(self, message=''):
        "Log an info message"
        if self.log_level <= dolfin.INFO:
            self._write(message)
    
    def debug(self, message=''):
        "Log a debug message"
        if self.log_level <= dolfin.DEBUG:
            self._write(message)
    
    def setup(self):
        """
        Setup logging to file if requested in the simulation input
        """
        log_name = self.simulation.input.get_output_file_path('output/log_name', None)
        if log_name is None:
            self.write_log = False
        else:
            self.write_log = True
            self.log_file_name = log_name
            self.log_file = open(self.log_file_name, 'wt')
    
    def end_of_simulation(self):
        """
        The simulation is done. Make sure the output
        file is flushed, but keep it open in case 
        some more output is coming
        """
        if self.write_log:
            self.log_file.flush()
