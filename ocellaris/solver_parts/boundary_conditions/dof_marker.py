import dolfin


def get_dof_region_marks(simulation, V):
    """
    Given a function space, return a dictionary mapping dof number to
    region number. Many dofs will not be included in the mapping since
    they are not inside a boundary region (not on a boundary facet)
    """
    pass
