import numpy
import matplotlib
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
import dolfin

class Plot2DDG0(object):
    def __init__(self, mesh, func):
        self.mesh = mesh
        self.func = func
        self.dofmap = func.function_space().dofmap().dofs()
        
        coords = mesh.coordinates()
        self.Ncell = mesh.num_entities(2)
        
        # Build matplotlib plolygon data
        self.verts = []
        for i, cell in enumerate(dolfin.cells(mesh)):
            corners = cell.entities(0)
            xs = [coords[i,0] for i in corners]
            ys = [coords[i,1] for i in corners]
            self.verts.append(zip(xs, ys))
    
    def plot(self, filename):
        values = numpy.zeros(self.Ncell, float)
        funcvals = self.func.vector()
        for i in xrange(self.Ncell):
            values[self.dofmap[i]] = funcvals[i]
        
        verts2, vals2 = [], []
        for vx, v in zip(self.verts, values):
            if v != 0:
                verts2.append(vx)
                vals2.append(v)
        vals2 = numpy.array(vals2) 
        
        # Make plot
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        polys = PolyCollection(verts2, array=vals2,
                               cmap=matplotlib.cm.jet,
                               edgecolors='black',
                               linewidths=0.25)
        ax.add_collection(polys)
        ax.autoscale_view()
        fig.colorbar(polys, ax=ax)
        
        fig.savefig(filename)
        pyplot.close(fig)
