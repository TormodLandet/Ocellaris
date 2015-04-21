import numpy
import dolfin
from .plot_CR1_2D_scalar import plot_2d_CR1

class PlotFacetExpressionDG0(object):
    def __init__(self, simulation, func, **options):
        """
        A plotter for FacetExpressionDG0 functions in 2D
        """
        self.func = func
        self.options = options
        
        # Get information about the underlying function space
        assert hasattr(func, 'ocellaris_cpp_expression_type')
        family = func.ocellaris_cpp_expression_type
        ndim = func.mesh.geometry().dim()
        
        # Check that the function is of a supported type
        assert family == 'FacetExpressionDG0'
        assert ndim == 2
        
        self.mesh = func.mesh
        
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
        scalars = self.func.facet_data.array() 
        radius = self.facet_length/8
        cmap = self.options.get('cmap', 'Reds')
        
        plot_2d_CR1(self.verts, self.facet_midpoints, radius, scalars, filename,
                    xlim=(self.coords[:,0].min(), self.coords[:,0].max()),
                    ylim=(self.coords[:,1].min(), self.coords[:,1].max()),
                    cmap=cmap) 
