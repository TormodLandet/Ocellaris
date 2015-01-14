import json
import collections
from .plot import Plotter
from .utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data

class Simulation(object):
    def __init__(self):
        self.input = {}
        self.data = {}
        self.plotting = Plotting(self)
        
        self.ndim = 0
        self.timestep = 0
        self.time = 0.0
        self.dt = 0.0
    
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
    
    def new_timestep(self, timestep_number, t, dt):
        """
        Several parts of the code wants to know these things,
        so we keep them in a central place
        """
        self.timestep = timestep_number
        self.time = t
        self.dt = dt
    
    def read_json_input_file(self, filename):
        with open(filename, 'rt') as inpf:
            inp = json.load(inpf, object_pairs_hook=collections.OrderedDict)
        assert inp['program'] == 'dgvof'
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
