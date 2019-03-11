# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from ..slope_limiter import SlopeLimiter


def create_component_limiters(simulation, vel_name, vel, comp_method):
    """
    Create component limiters for a velocity field. This ensures that
    the component input dictionary is populated with values from the
    parent velocity input unless explicitly overridden by the user
    """
    lims = []
    for d in range(simulation.ndim):
        comp_name = '%s%d' % (vel_name, d)

        vel_inp = simulation.input.ensure_path('slope_limiter/%s' % vel_name)
        comp_inp = simulation.input.ensure_path('slope_limiter/%s' % comp_name)

        # Make sure the default component limiter method is used if nothing
        # else is specified for this velocity component
        if 'method' not in comp_inp:
            comp_inp['method'] = comp_method

        # Copy (almost) all keys from the velocity to the component inputs
        # as long as they are not specifically given in the component input
        ignore = {'prelimiter', 'comp_method'}
        for key, value in vel_inp.items():
            if key not in comp_inp and key not in ignore:
                comp_inp[key] = value

        # Create the slope limiter for the component
        lim = SlopeLimiter(simulation, comp_name, vel[d])
        lims.append(lim)
    return lims
