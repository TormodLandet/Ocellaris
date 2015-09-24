import numpy
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
import dolfin
from ocellaris.utils import facet_dofmap

class Plot2DCR1Vec(object):
    def __init__(self, simulation, func, **options):
        """
        A plotter for Crouzeix-Raviart vector functions in 2D
        """
        self.func = func
        self.options = options
        
        # Get information about the underlying function space
        function_space = func.function_space()
        element = function_space.ufl_element()
        cell = function_space.ufl_cell()
        ndim = cell.topological_dimension()
        
        # Check that the function is of a supported type
        assert element.family() == 'Crouzeix-Raviart'
        assert element.degree() == 1
        assert ndim == 2
        assert function_space.num_sub_spaces() == 2
        
        self.mesh = function_space.mesh()
        self.dofmap0 = facet_dofmap(function_space.sub(0))
        self.dofmap1 = facet_dofmap(function_space.sub(1))
        
        # Build matplotlib plolygon data
        self.coords = self.mesh.coordinates()
        self.Ncell = self.mesh.num_entities(2)
        self.verts = []
        for cell in dolfin.cells(self.mesh):
            corners = cell.entities(0)
            xs = [self.coords[i,0] for i in corners]
            ys = [self.coords[i,1] for i in corners]
            self.verts.append(zip(xs, ys))
        
        # Build facet midpoint information
        self.Nfacet = self.mesh.num_entities(1)
        self.facet_centroids = []
        for facet in dolfin.facets(self.mesh):
            mp = facet.midpoint()
            self.facet_centroids.append((mp.x(), mp.y()))
        self.facet_centroids = numpy.array(self.facet_centroids, float)
    
    def plot(self, filename):
        """
        Plot the current state of the referenced function to a file
        """
        # Rearange function values according to the dofmap
        vectors = numpy.zeros((self.Nfacet, 2), float)
        funcvals = self.func.vector()
        for i, facet in enumerate(dolfin.facets(self.mesh)):
            fidx = facet.index()
            vectors[i, 0] = funcvals[self.dofmap0[fidx]]
            vectors[i, 1] = funcvals[self.dofmap1[fidx]]
        
        plot_2d_CR1_vec(self.verts, self.facet_centroids, vectors, filename,
                        xlim=(self.coords[:,0].min(), self.coords[:,0].max()),
                        ylim=(self.coords[:,1].min(), self.coords[:,1].max())) 
    
def plot_2d_CR1_vec(vertices, facet_centroids, vectors, filename, xlim, ylim):
    """
    Helper function to plot CR1 2D vector results. 
    Use the Plot2DCR1Vec class for a nicer interface to this function
    """
    # Make plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    
    polys = PolyCollection(vertices, facecolors='white', edgecolors='black', linewidths=0.25)
    ax.add_collection(polys)
    
    gx = vectors[:,0]
    gy = vectors[:,1]
    maxlen = (gx**2 + gy**2).max()**0.5
    target_maxlen = 0.1
    if maxlen == 0:
        maxlen = 1
    
    ax.quiver(facet_centroids[:,0], facet_centroids[:,1], gx+0.001, gy+0.001,
              angles='xy', scale_units='xy', scale=maxlen/target_maxlen)
    ax.axis('equal')
    
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    
    fig.tight_layout()
    fig.savefig(filename)
    pyplot.close(fig)
