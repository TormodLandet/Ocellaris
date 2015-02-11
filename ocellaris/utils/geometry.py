import numpy
import dolfin
from collections import namedtuple

# Compact way to store the information, may turn into classes later if needed
CellInfo = namedtuple('CellData', 'volume midpoint')
FacetInfo = namedtuple('FacetData', 'area midpoint normal')

def init_connectivity(simulation):
    """
    Initialize the needed connectivity data
    """
    mesh = simulation.data['mesh']
    ndim = simulation.ndim
    
    if ndim == 2:
        # Connectivity from face to edge
        mesh.init(2, 1)
        con21 = mesh.topology()(2, 1)
    
        # Connectivity from edge to face
        mesh.init(1, 2)
        con12 = mesh.topology()(1, 2)
        
        # Connectivity from face to face
        mesh.init(2, 2)
        con22 = mesh.topology()(2, 2)
        
        #simulation.data['connectivity_21'] = con21
        #simulation.data['connectivity_12'] = con12
        simulation.data['connectivity_FC'] = con12
        simulation.data['connectivity_CF'] = con21
        simulation.data['connectivity_CC'] = con22
    
    else:
        # Connectivity from cell to face
        mesh.init(3, 2)
        con32 = mesh.topology()(3, 2)
    
        # Connectivity from face to cell
        mesh.init(2, 3)
        con23 = mesh.topology()(2, 3)
        
        # Connectivity from cell to cell
        mesh.init(3, 3)
        con33 = mesh.topology()(3, 3)
        
        simulation.data['connectivity_FC'] = con23
        simulation.data['connectivity_CF'] = con32
        simulation.data['connectivity_CC'] = con33
    
def precompute_cell_data(simulation):
    """
    Get cell volume and midpoint in an easy to use format
    """
    mesh = simulation.data['mesh']
    ndim = simulation.ndim
    
    cell_info = {}        
    for cell in dolfin.cells(mesh):
        mp = cell.midpoint()
        if ndim == 2:
            midpoint = numpy.array([mp.x(), mp.y()], float)
        else:
            midpoint = numpy.array([mp.x(), mp.y(), mp.z()], float)
        
        volume = cell.volume()
        cell_info[cell.index()] = CellInfo(volume, midpoint)
    
    simulation.data['cell_info'] = cell_info

def precompute_facet_data(simulation):
    """
    Get facet normal and areas in an easy to use format
    """
    mesh = simulation.data['mesh']
    conFC = simulation.data['connectivity_FC']
    ndim = simulation.ndim
    cell_info = simulation.data['cell_info']
    
    # Get the facet areas from the cells
    areas = {}
    for cell in dolfin.cells(mesh):
        # Get the connected facets 
        if ndim == 2:
            facet_idxs = cell.entities(1)
        else:
            facet_idxs = cell.entities(2)
        
        # Loop over connected facets and get the area
        for i, fidx in enumerate(facet_idxs):
            a = cell.facet_area(i)
            if fidx in areas:
                assert a == areas[fidx]
            else:
                areas[fidx] = a
    
    # Loop over facets and gather the required information
    facet_info = {}
    for facet in dolfin.facets(mesh):
        fidx = facet.index()
        
        mp = facet.midpoint()
        if ndim == 2:
            midpoint = numpy.array([mp.x(), mp.y()], float)
        else:
            midpoint = numpy.array([mp.x(), mp.y(), mp.z()], float)
        
        # Find one cell connected to this facet. There can be one or two
        # connected cells, we only need the first one 
        connected_cells = conFC(fidx)
        icell0 = connected_cells[0]
        
        # Midpoint of local cell 0
        cell0_mp = cell_info[icell0].midpoint

        # Vector from local cell midpoint to face midpoint
        vec0 = midpoint - cell0_mp

        # Find a normal pointing out from cell 0
        normalpt = facet.normal()
        normal = numpy.array([normalpt.x(), normalpt.y()], float)
        if numpy.dot(vec0, normal) < 0:
            normal *= -1
        
        area = areas[fidx]
        facet_info[fidx] = FacetInfo(area, midpoint, normal)
    
    simulation.data['facet_info'] = facet_info
