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
    facet_dofmap = [None]*mesh.num_facets()
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

