import numpy
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
import dolfin

class Plot2DDG0Vec(object):
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
        
        # Check that the function is of a supported type
        assert element.family() == 'Discontinuous Lagrange'
        assert element.degree() == 0
        assert cell.topological_dimension() == 2
        assert function_space.num_sub_spaces() == 2
        
        self.mesh = function_space.mesh()
        self.dofmap0 = function_space.sub(0).dofmap().dofs()
        self.dofmap1 = function_space.sub(1).dofmap().dofs()
        
        # Build matplotlib polygon data
        self.coords = self.mesh.coordinates()
        self.Ncell = self.mesh.num_entities(2)
        self.verts = []
        self.centroids = numpy.zeros((self.Ncell, 2), float)
        for ic, cell in enumerate(dolfin.cells(self.mesh)):
            corners = cell.entities(0)
            xs = [self.coords[i,0] for i in corners]
            ys = [self.coords[i,1] for i in corners]
            self.verts.append(zip(xs, ys))
            c = cell.midpoint()
            self.centroids[ic] = (c.x(), c.y())
    
    def plot(self, filename, skip_zero_values=False):
        """
        Plot the current state of the referenced function to a file
        """
        # Rearange function values according to the dofmap
        vectors = numpy.zeros((self.Ncell, 2), float)
        funcvals = self.func.vector()
        for i, cell in enumerate(dolfin.cells(self.mesh)):
            cidx = cell.index()
            vectors[i,0] = funcvals[self.dofmap0[cidx]]
            vectors[i,1] = funcvals[self.dofmap1[cidx]]
        
        plot_2d_DG0_vec(self.verts, self.centroids, vectors, filename,
                        xlim=(self.coords[:,0].min(), self.coords[:,0].max()),
                        ylim=(self.coords[:,1].min(), self.coords[:,1].max())) 
    
def plot_2d_DG0_vec(vertices, centroids, vectors, filename, xlim, ylim):
    """
    Helper function to plot DG0 2D vector results. 
    Use the Plot2DDG0Vec class for a nicer interface to this function
    """
    # Make plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    polys = PolyCollection(vertices, facecolor='white', edgecolors='black', linewidths=0.25)
    
    ax.add_collection(polys)

    gx = vectors[:,0]
    gy = vectors[:,1]
    maxlen = (gx**2 + gy**2).max()**0.5
    target_maxlen = 0.1
    if maxlen == 0:
        maxlen = 1
    
    ax.quiver(centroids[:,0], centroids[:,1], gx+0.001, gy+0.001,
              angles='xy', scale_units='xy', scale=maxlen/target_maxlen)
    ax.axis('equal')
    
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    
    fig.tight_layout()
    fig.savefig(filename)
    pyplot.close(fig)
