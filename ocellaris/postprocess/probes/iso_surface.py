import dolfin
import numpy
from matplotlib import pyplot
from ocellaris.utils import gather_lines_on_root, timeit, ocellaris_error
from . import Probe, register_probe


@register_probe('IsoSurface')
class IsoSurface(Probe):
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.cells_with_surface = None
        inp = probe_input
        
        assert self.simulation.ndim == 2, 'IsoSurface only implemented in 2D (contour line)'

        # Read input
        self.name = inp.get_value('name', required_type='string')
        self.field_name = inp.get_value('field', required_type='string')
        self.value = inp.get_value('value', required_type='float')
        self.custom_hook_point = inp.get_value('custom_hook', None, required_type='string')
        self.field = simulation.data[self.field_name]
        
        # Should we write the data to a file
        prefix = simulation.input.get_value('output/prefix', None, 'string')
        file_name = inp.get_value('file_name', None, 'string')
        self.write_file = file_name is not None
        if self.write_file:
            if prefix is not None:
                self.file_name = prefix + file_name
            else:
                self.file_name = file_name
            self.write_interval = inp.get_value('write_interval', 1, 'int')
        
        # Should we pop up a matplotlib window when running?
        self.show_interval = inp.get_value('show_interval', 0, 'int')
        self.show = self.show_interval != 0 and simulation.rank == 0
        self.xlim = inp.get_value('xlim', (None, None), 'list(float)')
        self.ylim = inp.get_value('ylim', (None, None), 'list(float)')
        if not len(self.xlim) == 2:
            ocellaris_error('Plot xlim must be two numbers',
                            'IsoSurface probe "%s" contains invalid xlim specification' % self.name)
        if not len(self.ylim) == 2:
            ocellaris_error('Plot ylim must be two numbers',
                            'IsoSurface probe "%s" contains invalid ylim specification' % self.name)
        
        if self.write_file and simulation.rank == 0:
            self.output_file = open(self.file_name, 'wt')
            self.output_file.write('# Ocellaris iso surface of the %s field\n' % self.field_name)
            self.output_file.write('# value = %15.5e\n' % self.value)
            self.output_file.write('# dim = %d\n' % self.simulation.ndim)
        
        if self.show and simulation.rank == 0:
            pyplot.ion()
            self.fig = pyplot.figure()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_title('Iso surface %s' % self.name)
        
        if self.custom_hook_point is not None:
            simulation.hooks.add_custom_hook(self.custom_hook_point, self.run, 'Probe "%s"' % self.name)
        else:
            self.end_of_timestep = self.run
    
    def run(self):
        """
        Find and output the line probe
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
        surfaces, cells = get_iso_surfaces(self.simulation, self.field, self.value)
        self.cells_with_surface = cells
        
        # Create lines (this assumes 2D and throws away the z-component)
        lines = []
        for surface in surfaces:
            x = numpy.array([pos[0] for pos in surface], float)
            y = numpy.array([pos[1] for pos in surface], float)
            lines.append((x, y))
        
        # Communicate lines to the root process in case we are running in parallel
        gather_lines_on_root(lines)
        
        if update_file and self.simulation.rank == 0:
            self.output_file.write('Time %10.5f nsurf %d\n' % (self.simulation.time, len(lines)))
            for x, y in lines:
                self.output_file.write(' '.join('%10.5f' % v for v in x) + '\n')
                self.output_file.write(' '.join('%10.5f' % v for v in y) + '\n')
                self.output_file.write(' '.join('%10.5f' % 0 for v in x) + '\n')
        
        if update_plot and self.simulation.rank == 0:
            self.ax.clear()
            for x, y in lines:
                self.ax.plot(x, y)
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
        if self.write_file and self.simulation.rank == 0:
            self.output_file.close()


@timeit
def get_iso_surfaces(simulation, field, value):
    """
    Find the iso-surfaces (contour lines) of the
    given field with the given scalar value 
    """
    assert simulation.ndim == 2
    
    if field.function_space().ufl_element().degree() == 0:
        return get_iso_surfaces_picewice_constants(simulation, field, value)
    
    mesh = simulation.data['mesh']
    all_values = field.compute_vertex_values()
    
    # We will collect the cells containing the iso surface
    cells_with_surface = numpy.zeros(mesh.num_cells(), bool)
    connectivity_FC = simulation.data['connectivity_FC']
    
    # Find the crossing points where the contour crosses a facet
    crossing_points = {}
    for facet in dolfin.facets(mesh):
        fid = facet.index()
        
        # Get connected vertices and the field values there
        vertex_coords = []
        vertex_values = []
        for vertex in dolfin.vertices(facet):
            pt = vertex.point()
            vertex_coords.append((pt.x(), pt.y(), pt.z()))
            vertex_values.append(all_values[vertex.index()])
        assert len(vertex_coords) == 2
        
        # Check for iso surface crossing
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
        crossing_points[fid] = (x, y, z)
        
        # Find the cells connected to this facet
        for cid in connectivity_FC(fid):
            cells_with_surface[cid] = True
    
    # Get facet-facet connectivity via cells
    conFC = simulation.data['connectivity_FC']
    conCF = simulation.data['connectivity_CF']
    
    # Find facet to facet connections
    connections = {}
    for facet_id in crossing_points:
        connections[facet_id] = []
        for cell_id in conFC(facet_id):
            for facet_neighbour_id in conCF(cell_id):
                if facet_neighbour_id != facet_id and facet_neighbour_id in crossing_points:
                    connections[facet_id].append(facet_neighbour_id)
    
    # Make continous contour lines
    # Find end points of contour lines and start with these
    end_points = [facet_id for facet_id, neighbours in connections.items() if len(neighbours) == 1]
    contours_from_endpoints = contour_lines_from_endpoints(end_points, crossing_points, connections)
    
    # Include crossing points without neighbours or joined circles without end points
    other_points = crossing_points.keys()
    contours_from_singles_and_loops = contour_lines_from_endpoints(other_points, crossing_points, connections)
    
    assert len(crossing_points) == 0
    return contours_from_endpoints + contours_from_singles_and_loops, cells_with_surface


def get_iso_surfaces_picewice_constants(simulation, field, value):
    """
    Find the iso-surfaces (contour lines) of the
    given field with the given scalar value
    
    The field is assumed to be piecewice constant (DG0) 
    """
    mesh = simulation.data['mesh']
    all_values = field.vector().get_local()
    dofmap = field.function_space().dofmap()
    startdof = dofmap.ownership_range()[0]
    
    # Mesh connectivities
    conFC = simulation.data['connectivity_FC']
    conVF = simulation.data['connectivity_VF']
    conFV = simulation.data['connectivity_FV']
    
    # We will collect the cells containing the iso surface
    cells_with_surface = numpy.zeros(mesh.num_cells(), bool)
    
    # We define acronym LCCM: line connecting cell midpoints
    #   - We restrinct ourselves to LCCMs that cross only ONE facet
    #   - We number LLCMs by the index of the crossed facet
    
    # Find the crossing points where the contour crosses a LCCM
    vertex_coords = numpy.zeros((2, 3), float)
    vertex_values = numpy.zeros(2, float) 
    crossing_points = {}
    for facet in dolfin.facets(mesh):
        fid = facet.index()
        cell_ids = conFC(fid)
        if len(cell_ids) != 2:
            continue
        
        has_ghost_cell = False
        for i, cell_id in enumerate(cell_ids):
            cell = dolfin.Cell(mesh, cell_id)
            if cell.is_ghost():
                has_ghost_cell = True
                break
            
            # LCCM endpoint coordinates
            pt = cell.midpoint()
            vertex_coords[i,0] = pt.x()
            vertex_coords[i,1] = pt.y()
            vertex_coords[i,2] = pt.z()
            
            # LCCM endpoint values
            dofs = dofmap.cell_dofs(cell_id)
            assert len(dofs) == 1
            vertex_values[i] = all_values[dofs[0] - startdof]
            
        if has_ghost_cell:
            continue
        
        b1, b2 = vertex_values[0] < value, vertex_values[1] < value
        if (b1 and b2) or not (b1 or b2):
            # LCCM not crossed by contour
            continue
        
        # Find the location where the contour line crosses the LCCM
        v1, v2 = vertex_values
        fac = (v1 - value)/(v1 - v2)
        crossing_points[fid] = (1 - fac)*vertex_coords[0] + fac*vertex_coords[1]
        
        # Find the cell containing the contour line
        surf_cid = cell_ids[0] if fac <= 0.5 else cell_ids[1]
        cells_with_surface[surf_cid] = True
    
    # Find facet to facet connections
    connections = {}
    for facet_id in crossing_points:
        connections[facet_id] = []
        for vertex_id in conFV(facet_id):
            for facet_neighbour_id in conVF(vertex_id):
                if facet_neighbour_id != facet_id and facet_neighbour_id in crossing_points:
                    connections[facet_id].append(facet_neighbour_id)
    
    # Make continous contour lines
    # Find end points of contour lines and start with these
    end_points = [facet_id for facet_id, neighbours in connections.items() if len(neighbours) == 1]
    contours_from_endpoints = contour_lines_from_endpoints(end_points, crossing_points, connections)
    
    # Include crossing points without neighbours or joined circles without end points
    other_points = crossing_points.keys()
    contours_from_singles_and_loops = contour_lines_from_endpoints(other_points, crossing_points, connections)
    
    assert len(crossing_points) == 0
    return contours_from_endpoints + contours_from_singles_and_loops, cells_with_surface


def contour_lines_from_endpoints(endpoints, crossing_points, connections):
    """
    Given facet ids of endpoints, follow the contour line and create contours
    """
    contours = []
    for endpoint in endpoints:
        if not endpoint in crossing_points:
            # This has been taken by the other end
            continue
        
        # Make a new contour line
        contour = [crossing_points.pop(endpoint)]
        
        # Loop over neighbours to the end of the contour
        queue = list(connections[endpoint])
        prev = endpoint
        while queue:
            facet_id = queue.pop()
            if facet_id in crossing_points:
                contour.append(crossing_points.pop(facet_id))
                queue.extend(connections[facet_id])
                prev = facet_id
            if facet_id == endpoint and prev in connections[endpoint] and len(contour) != 2:
                contour.append(contour[0])
                break
        
        contours.append(contour)
    
    return contours
