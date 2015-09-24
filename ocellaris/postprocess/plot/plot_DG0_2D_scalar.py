import numpy
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
import dolfin

class Plot2DDG0(object):
    def __init__(self, simulation, func, **options):
        """
        A plotter for DG0 functions in 2D
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
        assert element.degree() == 0
        assert cell.topological_dimension() == 2
        
        self.mesh = function_space.mesh()
        self.dofmap = function_space.dofmap().dofs()
        
        # Build matplotlib polygon data
        self.coords = self.mesh.coordinates()
        self.Ncell = self.mesh.num_entities(2)
        self.verts = []
        for cell in dolfin.cells(self.mesh):
            corners = cell.entities(0)
            xs = [self.coords[i,0] for i in corners]
            ys = [self.coords[i,1] for i in corners]
            self.verts.append(zip(xs, ys))
    
    def plot(self, filename, skip_zero_values=False):
        """
        Plot the current state of the referenced function to a file
        """
        # Rearange function values according to the dofmap
        values = numpy.zeros(self.Ncell, float)
        funcvals = self.func.vector()
        for i in xrange(self.Ncell):
            values[self.dofmap[i]] = funcvals[i]
        
        # Only include polygons with non-zero values
        if skip_zero_values:
            verts2, vals2 = [], []
            for vx, v in zip(self.verts, values):
                if v != 0:
                    verts2.append(vx)
                    vals2.append(v)
            vals2 = numpy.array(vals2)
        else:
            verts2 = self.verts
            vals2 = values
        
        clim = self.options.get('clim', None)
        
        plot_2d_DG0(verts2, vals2, filename,
                    xlim=(self.coords[:,0].min(), self.coords[:,0].max()),
                    ylim=(self.coords[:,1].min(), self.coords[:,1].max()),
                    clim=clim) 
    
def plot_2d_DG0(vertices, values, filename, xlim, ylim, clim=None, cmap='Blues'):
    """
    Helper function to plot DG0 2D results. 
    Use the Plot2DDG0 class for a nicer interface to this function
    """
    # Make plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    polys = PolyCollection(vertices, array=values,
                           cmap=cmap,
                           edgecolors='black',
                           linewidths=0.25)
    
    if clim is not None:
        polys.set_clim(*clim)
        polys.cmap.set_under('grey')
        polys.cmap.set_over('grey')
    
    ax.add_collection(polys)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fig.colorbar(polys, ax=ax)
    
    fig.tight_layout()
    fig.savefig(filename)
    pyplot.close(fig)
