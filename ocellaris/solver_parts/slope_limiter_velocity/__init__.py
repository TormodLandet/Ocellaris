# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from ocellaris.utils import ocellaris_error, verify_key


LIMITER = 'None'
USE_CPP = True
_VELOCITY_SLOPE_LIMITERS = {}


def add_velocity_slope_limiter(name, slope_limiter_class):
    """
    Register a slope limiter
    """
    _VELOCITY_SLOPE_LIMITERS[name] = slope_limiter_class


def register_velocity_slope_limiter(name):
    """
    A class decorator to register slope limiters
    """

    def register(slope_limiter_class):
        add_velocity_slope_limiter(name, slope_limiter_class)
        return slope_limiter_class

    return register


def get_velocity_slope_limiter(name):
    """
    Return a slope limiter by name
    """
    try:
        return _VELOCITY_SLOPE_LIMITERS[name]
    except KeyError:
        ocellaris_error(
            'Velocity slope limiter "%s" not found' % name,
            'Available velocity slope limiters:\n'
            + '\n'.join(
                '  %-20s - %s' % (n, s.description)
                for n, s in sorted(_VELOCITY_SLOPE_LIMITERS.items())
            ),
        )
        raise


class VelocitySlopeLimiterBase(object):
    description = 'No description available'
    additional_plot_funcs = ()
    active = True


@register_velocity_slope_limiter('None')
class DoNothingSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'No limiting'
    active = False

    def __init__(self, *argv, **kwargs):
        pass

    def run(self):
        pass


def SlopeLimiterVelocity(
    simulation,
    vel_u,
    vel_name,
    default_limiter=LIMITER,
    vel_w=None,
    default_use_cpp=USE_CPP,
):
    """
    Limit the slope of the given vector field to obtain boundedness
    """
    # Get user provided input (or default values)
    inp = simulation.input.get_value('slope_limiter/%s' % vel_name, {}, 'Input')
    method = inp.get_value('method', default_limiter, 'string')
    use_cpp = inp.get_value('use_cpp', default_use_cpp, 'bool')
    plot_exceedance = inp.get_value('plot', False, 'bool')

    # Get the limiter
    simulation.log.info(
        '    Using velocity slope limiter %s for %s' % (method, vel_name)
    )
    limiter_class = get_velocity_slope_limiter(method)
    limiter = limiter_class(simulation, vel_u, vel_name, vel_w, use_cpp)

    # Add extra limiter outputs
    if plot_exceedance:
        for func in limiter.additional_plot_funcs:
            simulation.io.add_extra_output_function(func)

    return limiter


from . import solenoidal_limiter
from . import componentwise_limiter
