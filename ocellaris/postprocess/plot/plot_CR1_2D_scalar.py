from math import pi
import numpy
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
import dolfin
from ocellaris.utils import facet_dofmap

class Plot2DCR1(object):
    def __init__(self, simulation, func, **options):
        """
        A plotter for Crouzeix-Raviart functions in 2D
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
        assert function_space.num_sub_spaces() == 0
        
        self.mesh = function_space.mesh()
        self.dofmap = facet_dofmap(function_space)
        
        # Build matplotlib plolygon data
        self.coords = self.mesh.coordinates()
        self.Ncell = self.mesh.num_entities(2)
        self.verts = []
        for cell in dolfin.cells(self.mesh):
            corners = cell.entities(0)
            xs = [self.coords[i,0] for i in corners]
            ys = [self.coords[i,1] for i in corners]
            self.verts.append(zip(xs, ys))
        
        # Create list of facet midpoints
        # and find the average facet length
        self.Nfacet = self.mesh.num_entities(1)
        facet_info = simulation.data['facet_info']
        self.facet_midpoints = numpy.zeros((self.Nfacet, 2), float)
        self.facet_length = 0
        for i, facet in enumerate(dolfin.facets(self.mesh)):
            fidx = facet.index()
            info = facet_info[fidx]
            self.facet_midpoints[i] = info.midpoint
            self.facet_length += info.area
        self.facet_length /= self.Nfacet
    
    def plot(self, filename):
        """
        Plot the current state of the referenced function to a file
        """
        # Rearange function values according to the dofmap
        scalars = numpy.zeros(self.Nfacet, float)
        funcvals = self.func.vector()
        
        for i, facet in enumerate(dolfin.facets(self.mesh)):
            fidx = facet.index()
            scalars[i] = funcvals[self.dofmap[fidx]]
        
        radius = self.facet_length/8
        cmap = self.options.get('cmap', 'Reds')
        
        plot_2d_CR1(self.verts, self.facet_midpoints, radius, scalars, filename,
                    xlim=(self.coords[:,0].min(), self.coords[:,0].max()),
                    ylim=(self.coords[:,1].min(), self.coords[:,1].max()),
                    cmap=cmap) 
    
def plot_2d_CR1(vertices, facet_centroids, radius, scalars, filename, xlim, ylim, cmap):
    """
    Helper function to plot CR1 2D scalar results. 
    Use the Plot2DCR1 class for a nicer interface to this function
    """
    # Make plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    
    def fix_axes():
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        fig.tight_layout()
    
    # Find  the length of the scalar circle radius in pixels
    fix_axes()
    ptx, pty = ax.transData.transform([(radius, 0), (0, radius)]) - ax.transData.transform((0,0))
    rad_px = min(ptx[0], pty[1])
    area_px = pi*rad_px**2
    
    vmin = scalars.min()
    vmax = scalars.max()
    if vmin == vmax == 0:
        vmax = 1
    
    # Plot the scalars with colors
    patches = ax.scatter(facet_centroids[:,0], facet_centroids[:,1], s=area_px, c=scalars, cmap=cmap,
                         vmin=vmin, vmax=vmax, linewidths=0)
    fig.colorbar(patches, ax=ax)
    fix_axes()
    
    # Add grid data above the scalars
    polys = PolyCollection(vertices, facecolors=(1, 1, 1, 0), edgecolors='black', linewidths=0.25)
    ax.add_collection(polys)
    
    fig.savefig(filename)
    pyplot.close(fig)
