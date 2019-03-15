# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import shift_fields


def before_simulation(simulation, force_steady=False):
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
    maxabs = dolfin.MPI.max(dolfin.MPI.comm_world, float(maxabs))
    if maxabs > 0:
        starting_order = 2

    if force_steady:
        simulation.log.info('Setting time derivatives to zero')
        simulation.data['time_coeffs'].assign(dolfin.Constant([0.0, 0.0, 0.0]))
    elif starting_order == 2:
        # Switch to second order time stepping
        simulation.log.info(
            'Initial values for upp are found and used, '
            'starting with second order time stepping.'
        )
        simulation.data['time_coeffs'].assign(dolfin.Constant([3 / 2, -2.0, 1 / 2]))
    else:
        # Standard first order time stepping
        simulation.log.info(
            'Initial values for upp are not found, ' 'starting with first order time stepping.'
        )
        simulation.data['time_coeffs'].assign(dolfin.Constant([1.0, -1.0, 0.0]))
    update_convection(simulation, starting_order, force_steady=force_steady)

    simulation.log.info('\nTime loop is now starting\n', flush='force')


def update_timestep(simulation):
    """
    Switch to first order time stepping if the timestep has changed
    """
    dt = simulation.input.get_value('time/dt', required_type='float')
    dt_prev = simulation.dt

    if dt != dt_prev:
        simulation.log.info('Temporarily changing to first order time stepping')
        simulation.data['time_coeffs'].assign(dolfin.Constant([1.0, -1.0, 0.0]))

    return dt


def after_timestep(simulation, is_steady, force_steady=False):
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

    if force_steady:
        simulation.data['time_coeffs'].assign(dolfin.Constant([0.0, 0.0, 0.0]))
    else:
        # Change time coefficient to second order
        simulation.data['time_coeffs'].assign(dolfin.Constant([3 / 2, -2, 1 / 2]))

    # Extrapolate the convecting velocity to the next step
    update_convection(simulation, force_steady=force_steady)

    return vel_diff


def update_convection(simulation, order=2, force_steady=False):
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

        if order == 1 or force_steady:
            uic.assign(uip)
        else:
            # Backwards difference formulation - standard linear extrapolation
            uic.vector().zero()
            uic.vector().axpy(2.0, uip.vector())
            uic.vector().axpy(-1.0, uipp.vector())
            uic.vector().apply('insert')
