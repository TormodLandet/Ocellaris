import dolfin

def create_krylov_solver(solver_name, preconditioner_name, parameters):
    """
    Create a Krylov solver with the given solver and preconditioner
    
    The parameters argument is a *list* of dictionaries which are
    to be used as parameters to the Krylov solver. Settings in the
    first dictionary in this list will be (potentially) overwritten
    by settings in later dictionaries. The use case is to provide
    sane defaults as well as allow the user to override the defaults
    in the input file  
    """
    solver = dolfin.KrylovSolver(solver_name, preconditioner_name)
    
    for parameter_set in parameters:
        apply_settings(solver.parameters, parameter_set)
        
    return solver

def apply_settings(parameters, new_values):
    """
    This functiuon does almost the same as::
    
        parameters.update(new_values)
        
    The difference is that subdictionaries are handled
    recursively and not replaced outright
    """
    for key, value in new_values.items():
        if isinstance(value, dict):
            apply_settings(parameters[key], value)
        else:
            parameters[key] = value
