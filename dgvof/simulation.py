import json
import collections
from .plot import Plotter
from .utils.geometry import init_connectivity, precompute_cell_data, precompute_facet_data

class Simulation(object):
    def __init__(self):
        self.input = {}
        self.data = {}
        self.plots = {}
    
    def set_mesh(self, mesh):
        self.data['mesh'] = mesh
        init_connectivity(self)
        precompute_cell_data(self)
        precompute_facet_data(self)
    
    def add_plot(self, plot_name, plotter):
        """
        Add a plot to the simulation
        """
        if not hasattr(plotter, 'plot'):
            # This is not a plotter but something that can be plotted
            plotter = Plotter(plotter)
        
        self.plots[plot_name] = plotter
    
    def plot_all(self, timestep, t):
        for name in self.plots:
            self.plot(name, timestep, t)
    
    def plot(self, name, timestep, t):
        name_template_dafault = 'fig/{name}_{timestep:07d}_{t:010.6f}.png'
        name_template = self.input.get('plots', {}).get('name_template', name_template_dafault)
        filename = name_template.format(name=name, timestep=timestep, t=t)
        self.plots[name].plot(filename) 
    
    def read_json_input_file(self, filename):
        with open(filename, 'rt') as inpf:
            inp = json.load(inpf, object_pairs_hook=collections.OrderedDict)
        assert inp['program'] == 'dgvof'
        assert inp['version'] == 1.0
        assert inp['type'] == 'input'
        self.input = inp

