from ocellaris.postprocess import Plotter


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
