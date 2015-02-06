from . import register_boundary_condition, BoundaryCondition

@register_boundary_condition('ConstantGradient')
class NeumannBC(BoundaryCondition):
    description = 'A prescribed constant value Neumann condition'
    
    def __init__(self, simulation, var_name, inp_dict, marker, mark_id):
        """
        Neumann condition
        """
        value = inp_dict['value']
        if isinstance(value, list):
            for d in range(simulation.ndim):
                pass
        else:
            pass
             
