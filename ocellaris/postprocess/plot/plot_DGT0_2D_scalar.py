import numpy
import dolfin
from ocellaris.utils import facet_dofmap
from .plot_CR1_2D_scalar import plot_2d_CR1

class PlotDGT0(object):
    def __init__(self, simulation, func, **options):
        """
        A plotter for DGT0 functions in 2D
        """
        self.func = func
        self.options = options
        
        # Get information about the underlying function space
        function_space = func.function_space()
        family = func.ufl_element().family()
        mesh = function_space.mesh()
        ndim = mesh.geometry().dim()
        
        # Check that the function is of a supported type
        assert family == 'Discontinuous Lagrange Trace'
        assert ndim == 2
        
        self.mesh = mesh
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
