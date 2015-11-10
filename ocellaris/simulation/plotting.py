import os
from ocellaris.postprocess import Plotter


class Plotting(object):
    def __init__(self, simulation):
        """
        Central place to register plots that will be output during
        the simulation
        """
        self.simulation = simulation
        self.plots = {}
        self.active = simulation.rank == 0
    
    def add_plot(self, plot_name, plotter, **options):
        """
        Add a plot to the simulation
        """
        if not self.active:
            return
        
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
        if not self.active:
            return
        sim = self.simulation
        
        # Get the figure directory
        output_prefix = sim.input.get_value('output/prefix', '', 'string')
        figdir = output_prefix + '_plots'
        if not os.path.isdir(figdir):
            os.mkdir(figdir)
        
        # Get the file name
        name_template_dafault = '{figdir}/{name}_{timestep:07d}_{t:010.6f}{extra}.png'
        name_template = sim.input.get_value('plots/name_template', name_template_dafault, 'string')
        filename = name_template.format(figdir=figdir, name=name, timestep=sim.timestep, t=sim.time, extra=extra)
        
        # Create the image
        self.plots[name].plot(filename)
