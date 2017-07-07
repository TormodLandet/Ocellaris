import dolfin
from ocellaris.utils import shift_fields


def before_simulation(simulation):
    """
    Handle timestepping issues before starting the simulation. There are
    basically two options, either we have full velocity history available,
    either from initial conditions on the input file or from a restart file,
    or there is only access to one previous time step and we need to start
    up using first order timestepping
    """
    starting_order = 1
    
    # Check if there are non-zero values in the upp vectors
    maxabs = 0
    for d in range(simulation.ndim):
        this_maxabs = abs(simulation.data['upp%d' % d].vector().get_local()).max()
        maxabs = max(maxabs, this_maxabs)
    maxabs = dolfin.MPI.max(dolfin.mpi_comm_world(), float(maxabs))
    if maxabs > 0:
        starting_order = 2
    
    # Switch to second order time stepping
    if starting_order == 2:
        simulation.log.info('Initial values for upp are found and used')
        simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
    update_convection(simulation, starting_order)


def after_timestep(simulation, is_steady):
    """
    Move u -> up, up -> upp and prepare for the next time step
    """
    # Stopping criteria for steady state simulations  
    vel_diff = None
    if is_steady:
        vel_diff = 0
        for d in range(simulation.ndim):
            u_new = simulation.data['u%d' % d]
            up = simulation.data['up%d' % d]
            diff = abs(u_new.vector().get_local() - up.vector().get_local()).max()
            vel_diff = max(vel_diff, diff)
    
    shift_fields(simulation, ['u%d', 'up%d', 'upp%d'])
    shift_fields(simulation, ['u_conv%d', 'up_conv%d', 'upp_conv%d'])
    
    # Change time coefficient to second order
    simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
    
    # Extrapolate the convecting velocity to the next step
    update_convection(simulation)
    
    return vel_diff


def update_convection(simulation, order=2):
    """
    Update terms used to linearise and discretise the convective term
    """
    ndim = simulation.ndim
    data = simulation.data
    
    # Update convective velocity field components
    for d in range(ndim):
        uic = data['u_conv%d' % d]
        uip = data['up_conv%d' % d]
        uipp = data['upp_conv%d' % d]
        
        if order == 1:
            uic.assign(uip)
        else:
            # Backwards difference formulation - standard linear extrapolation
            uic.vector().zero()
            uic.vector().axpy(2.0, uip.vector())
            uic.vector().axpy(-1.0, uipp.vector())
            uic.vector().apply('insert')
