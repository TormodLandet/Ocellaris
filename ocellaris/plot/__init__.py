from .plot_DG0_2D_scalar import Plot2DDG0
from .plot_DG0_2D_vector import Plot2DDG0Vec
from .plot_CR1_2D_scalar import Plot2DCR1
from .plot_CR1_2D_vector import Plot2DCR1Vec

def Plotter(simulation, dolfin_function, *args, **kwargs):
    """
    Factory function for plotting
    """
    # Get information about the underlying function space
    function_space = dolfin_function.function_space()
    element = function_space.ufl_element()
    cell = function_space.cell()
    
    # Dispatch to the correct implementation class
    # TODO: make a registry here or ask the plotting classes if they support the function?
    function_type = (element.family(), element.degree(), cell.topological_dimension(), function_space.num_sub_spaces()) 
    
    if function_type == ('Discontinuous Lagrange', 0, 2, 0):
        return Plot2DDG0(simulation, dolfin_function, *args, **kwargs)
    if function_type == ('Discontinuous Lagrange', 0, 2, 2):
        return Plot2DDG0Vec(simulation, dolfin_function, *args, **kwargs)
    elif function_type == ('Crouzeix-Raviart', 1, 2, 0):
        return Plot2DCR1(simulation, dolfin_function, *args, **kwargs)
    elif function_type == ('Crouzeix-Raviart', 1, 2, 2):
        return Plot2DCR1Vec(simulation, dolfin_function, *args, **kwargs)
    else:
        raise ValueError('No plotters for "%s" of order %d in %dD with %d sub spaces' % function_type)
