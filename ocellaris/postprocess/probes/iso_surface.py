import dolfin
import numpy
from matplotlib import pyplot
from . import Probe, register_probe

@register_probe('IsoSurface')
class IsoSurface(Probe):
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.input = probe_input
        
        assert self.simulation.ndim == 2, 'IsoSurface only implemented in 2D (contour line)'

        # Read input
        name = self.input['name']
        self.field_name = self.input['field']
        self.value = self.input['value']
        
        self.field = simulation.data[self.field_name]
        
        # Should we write the data to a file
        prefix = simulation.input.get_value('output/prefix', None, 'string')
        file_name = self.input.get('file_name', '')
        self.write_file = file_name is not None
        if self.write_file:
            if prefix is not None:
                self.file_name = prefix + file_name
            else:
                self.file_name = file_name
            self.write_interval = self.input.get('write_interval', 1)
        
        # Should we pop up a matplotlib window when running?
        self.show_interval = self.input.get('show_interval', 0)
        self.show = self.show_interval != 0
        self.xlim = self.input.get('xlim', (None, None))
        self.ylim = self.input.get('ylim', (None, None))
        
        if self.write_file:
            self.output_file = open(self.file_name, 'wt')
            self.output_file.write('# Ocellaris iso surface of the %s field\n' % self.field_name)
            self.output_file.write('# value = %15.5e\n' % self.value)
            self.output_file.write('# dim = %d\n' % self.simulation.ndim)
        
        if self.show:
            pyplot.ion()
            self.fig = pyplot.figure()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_title('Iso surface %s' % name)
            self.line, = self.ax.plot([], [], 'ko')
    
    def end_of_timestep(self):
        """
        Output the line probe at the end of the
        """
        it = self.simulation.timestep
        
        # Should we update the plot?
        update_plot = False 
        if self.show and (it == 1 or it % self.show_interval == 0):
            update_plot = True

        # Should we update the file?
        update_file = False
        if self.write_file and (it == 1 or it % self.write_interval == 0):
            update_file = True
        
        # Do not do any postprocessing for non-requested time steps
        if not (update_file or update_plot):
            return
        
        # Get the iso surfaces
        surfaces = get_iso_surfaces(self.simulation, self.field, self.value)
        
        if update_file:
            self.output_file.write('Time %10.5f nsurf %d\n' % (self.simulation.time, len(surfaces)))
            for surface in surfaces:
                self.output_file.write(' '.join('%10.5f' % pos[0] for pos in surface) + '\n')
                self.output_file.write(' '.join('%10.5f' % pos[1] for pos in surface) + '\n')
                self.output_file.write(' '.join('%10.5f' % pos[2] for pos in surface) + '\n')
        
        if update_plot:
            xvals, yvals = [], []
            for surface in surfaces:
                xvals.extend(pos[0] for pos in surface)
                yvals.extend(pos[1] for pos in surface)
            self.line.set_xdata(xvals)
            self.line.set_ydata(yvals)
            self.ax.set_xlabel('x')
            self.ax.set_ylabel('y')
            self.ax.relim()
            self.ax.autoscale_view()
            if self.xlim != (None, None):
                self.ax.set_xlim(*self.xlim)
            if self.ylim != (None, None):
                self.ax.set_ylim(*self.ylim)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
    
    def end_of_simulation(self):
        """
        The simulation is done. Close the output file
        """
        if self.write_file:
            self.output_file.close()
            
def get_iso_surfaces(simulation, field, value):
    """
    Find the iso-surfaces (contour lines) of the
    given field with the given scalar value 
    """
    assert simulation.ndim == 2
    mesh = simulation.data['mesh']
    all_values = field.compute_vertex_values()
    
    # Find the crossing points where the contour crosses a facet
    crossing_points = []
    for facet in dolfin.facets(mesh):
        vertex_coords = []
        vertex_values = []
        for vertex in dolfin.vertices(facet):
            pt = vertex.point()
            vertex_coords.append((pt.x(), pt.y(), pt.z()))
            vertex_values.append(all_values[vertex.index()])
        assert len(vertex_coords) == 2
        
        b1, b2 = vertex_values[0] < value, vertex_values[1] < value
        if (b1 and b2) or not (b1 or b2):
            # Facet not crossed by contour
            continue
        
        # Find the location where the contour line crosses the facet
        v1, v2 = vertex_values
        fac = (v1 - value)/(v1 - v2)
        x = (1 - fac)*vertex_coords[0][0] + fac*vertex_coords[1][0]
        y = (1 - fac)*vertex_coords[0][1] + fac*vertex_coords[1][1]
        z = (1 - fac)*vertex_coords[0][2] + fac*vertex_coords[1][2]
        crossing_points.append((facet.index(), x, y, z))
        
    # Create continous contour lines
    contours = []
    for fidx, x, y, z in crossing_points:
        # Just make a new contour line per point for now
        contours.append([(x, y, z)])
    
    return contours
