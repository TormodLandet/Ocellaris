from __future__ import division
import numpy
import dolfin


def facet_dofmap(V):
    """
    When working with Crouzeix-Raviart and DGT elements with dofs on the facets
    it can be useful to get the dof corresponding to a facet index.
    This function returns a list which gives such a mapping.
    """
    mesh = V.mesh()
    dofmap = V.dofmap()
    
    ndim = V.ufl_cell().topological_dimension()
    
    # Loop through cells and get dofs for each cell
    facet_dofmap = [None] * mesh.num_facets()
    for cell in dolfin.cells(mesh):
        dofs = dofmap.cell_dofs(cell.index())
        
        if ndim == 2:
            facet_idxs = cell.entities(1)
        elif ndim == 3:
            facet_idxs = cell.entities(2)
        
        assert len(dofs) == len(facet_idxs)
        
        # Loop through connected facets and store dofs for each facet
        for fidx, dof in zip(facet_idxs, dofs):
            facet_dofmap[fidx] = dof
    
    return facet_dofmap


def get_dof_neighbours(V):
    """
    Given a DG function space find, for each dof, the indices
    of the cells with dofs at the same locations
    """
    dm = V.dofmap()
    gdim = V.mesh().geometry().dim()
    num_cells_all = V.mesh().num_cells()
    dof_coordinates = V.tabulate_dof_coordinates().reshape((-1, gdim))
    
    # Get "owning cell" indices for all dofs
    cell_for_dof = [None] * V.dim()
    for ic in xrange(num_cells_all):
        dofs = dm.cell_dofs(ic)
        for dof in dofs:
            assert cell_for_dof[dof] is None
            cell_for_dof[dof] = ic
    
    # Map dof coordinate to dofs, this is for DG so multiple dofs
    # will share the same location
    coord_to_dofs = {}
    max_neighbours = 0
    for dof in xrange(len(dof_coordinates)):
        coord = tuple(round(x, 5) for x in dof_coordinates[dof])
        dofs = coord_to_dofs.setdefault(coord, [])
        dofs.append(dof)
        max_neighbours = max(max_neighbours, len(dofs)-1)
    
    # Find number of neighbour cells and their indices for each dof
    num_neighbours = numpy.zeros(V.dim(), numpy.intc)
    neighbours = numpy.zeros((V.dim(), max_neighbours), numpy.intc) - 1
    for nbs in coord_to_dofs.values():
        # Loop through dofs at this location
        for dof in nbs:
            # Loop through the dofs neighbours
            for nb in nbs:
                # Skip the dof itself
                if dof == nb:
                    continue
                # Get the neighbours "owning cell" index and store this
                nb_cell = cell_for_dof[nb]
                nn_prev = num_neighbours[dof]
                neighbours[dof,nn_prev] = nb_cell
                num_neighbours[dof] += 1
    
    return num_neighbours, neighbours
