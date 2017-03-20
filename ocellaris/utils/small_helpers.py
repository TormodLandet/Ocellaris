import dolfin


def create_vector_functions(simulation, vec_name, comp_name, V):
    """
    Create vector functions and add them to the simulation
    Example::
    
        create_vector_functions(simulation, 'u', 'u%d', Vu)
    
    This will create simulation.data['u'] and simulation.data['uX']  
    """
    vec = []
    for d in range(simulation.ndim):
        name = comp_name % d
        f = dolfin.Function(V)
        f.rename(name, name)
        vec.append(f)
        simulation.data[name] = f
    simulation.data[vec_name] = dolfin.as_vector(vec)
    return simulation.data[vec_name]


def shift_fields(simulation, names):
    """
    Shift variables to next time step. Examples::
    
        shift_fields(sim, ['u%d', 'up%d', 'upp%d'])
        shift_fields(sim, ['rho', 'rho_p', 'rho_pp'])
    
    Will cause eg rho_pp = rho_p, rho_p = rho. This function can
    also be used for assigment, since a = b can be written::
    
        shift_fields(sim, ['b', 'a'])
    
    """
    # Detect whether we are treating a vector field or a scalar 
    ndim = None
    if '%d' in names[0]:
        ndim = simulation.ndim
    
    for i in range(len(names) - 1):
        target = names[-i-1]
        source = names[-i-2]
        if ndim is None:
            simulation.data[target].assign(simulation.data[source])
        else:
            for d in range(ndim):
                simulation.data[target % d].assign(simulation.data[source % d])
