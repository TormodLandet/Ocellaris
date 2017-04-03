import numpy


class SlopeLimiterBoundaryConditions(object):
    def __init__(self, simulation, field_name, boundary_marks, dim):
        """
        This class helps slope limiting of cells adjacent to the boundary by
        providing either values of the function and its derivatives at the
        boundary 
        """
        self.simulation = simulation
        self.field_name = field_name
        self.boundary_marks = boundary_marks
        self.dim = dim
        self.active = False
    
    def activate(self):
        self.active = True
    
    def get_bcs(self):
        boundary_dof_type = numpy.zeros(self.dim, numpy.intc)
        boundary_dof_value = numpy.zeros(self.dim, float)
        if not self.active:
            return boundary_dof_type, boundary_dof_value
        
        regions = self.simulation.data['boundary']
        
        return boundary_dof_type, boundary_dof_value
    