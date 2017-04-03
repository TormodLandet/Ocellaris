import numpy
import dolfin
from ocellaris.utils import timeit


@timeit
def get_dof_region_marks(simulation, V):
    """
    Given a function space, return a dictionary mapping dof number to
    region number. Many dofs will not be included in the mapping since
    they are not inside a boundary region (not on a boundary facet)
    """
    # This function only supports a small subset of function spaces
    family = V.ufl_element().family()
    degree = V.ufl_element().degree()
    assert family in ('Lagrange', 'Discontinuous Lagrange')
    assert V.mesh().topology().dim() == 2
    
    # Get local indices for the facet dofs for each facet in the cell
    if degree == 1:
        facet_dof_indices = numpy.zeros((3, 2), dtype=int)
        facet_dof_indices[0,:] = (1, 2)
        facet_dof_indices[1,:] = (0, 2)
        facet_dof_indices[2,:] = (0, 1)
    else:
        assert degree == 2
        facet_dof_indices = numpy.zeros((3, 3), dtype=int)
        facet_dof_indices[0,:] = (1, 2, 3)
        facet_dof_indices[1,:] = (0, 2, 4)
        facet_dof_indices[2,:] = (0, 1, 5)
    
    # Loop over mesh and get dofs that are connected to boundary regions
    dm = V.dofmap()
    facet_marks = [int(m) for m in simulation.data['boundary_marker'].array()]
    mesh = simulation.data['mesh']
    dof_region_marks = {}
    for cell in dolfin.cells(mesh):
        dofs = dm.cell_dofs(cell.index())
        
        for ifacet, facet in enumerate(dolfin.facets(cell)):
            # Get facet region marker
            mark = facet_marks[facet.index()] - 1
            
            # Skip non-boundary facets
            if mark == -1:
                continue
            
            facet_dofs = dofs[facet_dof_indices[ifacet]]
            for fdof in facet_dofs:
                dof_region_marks.setdefault(fdof, []).append(mark)
    
    assert len(dof_region_marks) > 4
    return dof_region_marks


def mark_cell_layers(simulation, V, layers=1, dof_region_marks=None):
    """
    Mark all dofs as boundary dofs if they are connected to a
    cell that contains boundary dofs. The initial dof list is
    taken as the keys from the given dof_region_mark dictionary  
    """
    if dof_region_marks is None:
        dof_region_marks = get_dof_region_marks(simulation, V)
    
    # Initialize the boolean dof marker
    marker = numpy.zeros(V.dim(), bool)
    for dof in dof_region_marks:
        marker[dof] =  True
    
    # Mark boundary cells
    mesh = simulation.data['mesh']
    dm = V.dofmap()    
    boundary_cells = set()
    for cell in dolfin.cells(mesh):
        dofs = dm.cell_dofs(cell.index())
        for dof in dofs:
            if marker[dof]:
                boundary_cells.add(cell.index())
                continue
    
    # Iteratively mark cells adjacent to boundary cells
    for _ in range(layers):
        boundary_cells_old = set(boundary_cells)
        for cell_index in boundary_cells_old:
            vertices = simulation.data['connectivity_CV'](cell_index)
            for vert_index in vertices:
                for nb in simulation.data['connectivity_VC'](vert_index):
                    boundary_cells.add(nb)
    
    for cell_index in boundary_cells:
        dofs = dm.cell_dofs(cell_index)
        for dof in dofs:
            marker[dof] = True
    
    return marker
