from matplotlib import pyplot
from matplotlib.tri import Triangulation

class Plot2DDG1(object):
    def __init__(self, simulation, func, **options):
        """
        A plotter for DG1 functions in 2D
        """
        self.func = func
        self.options = options
        
        # Get information about the underlying function space
        function_space = func.function_space()
        element = function_space.ufl_element()
        cell = function_space.ufl_cell()
        assert function_space.num_sub_spaces() == 0
        
        # Check that the function is of a supported type
        assert element.family() == 'Discontinuous Lagrange'
        assert element.degree() == 1
        assert cell.topological_dimension() == 2
        
        # Build matplotlib triangulation
        self.mesh = function_space.mesh()
        self.coords = self.mesh.coordinates()
        self.triangulation = Triangulation(self.coords[:, 0], self.coords[:, 1], self.mesh.cells())
        
    def plot(self, filename, skip_zero_values=False):
        """
        Plot the current state of the referenced function to a file
        """
        # TODO: this "projects" to CG1
        vertex_values = self.func.compute_vertex_values(self.mesh)
        clim = self.options.get('clim', None)
        
        plot_2d_DG1(self.triangulation, vertex_values, filename,
                    xlim=(self.coords[:,0].min(), self.coords[:,0].max()),
                    ylim=(self.coords[:,1].min(), self.coords[:,1].max()),
                    clim=clim) 
    
def plot_2d_DG1(triangulation, values, filename, xlim, ylim, clim=None, cmap='Blues'):
    """
    Helper function to plot DG1 2D results. 
    Use the Plot2DDG1 class for a nicer interface to this function
    """
    # Make plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    p = pyplot.tripcolor(triangulation, values, shading='gouraud')
    
    if clim is not None:
        p.set_clim(*clim)
        p.cmap.set_under('grey')
        p.cmap.set_over('grey')
    
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.colorbar(p, ax=ax)
    
    fig.tight_layout()
    fig.savefig(filename)
    pyplot.close(fig)
