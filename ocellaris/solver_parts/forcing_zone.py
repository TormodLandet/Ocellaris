# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import verify_field_variable_definition, verify_key


def add_forcing_zone(simulation, fzones, inp):
    """
    Add a penalty forcing zone to the simulation
    """
    name = inp.get_value('name', required_type='string')
    ztype = inp.get_value('type', required_type='string')
    penalty = inp.get_value('penalty', required_type='float')
    zone_vardef = inp.get_value('zone', required_type='string')
    target_vardef = inp.get_value('target', required_type='string')
    plot = inp.get_value('plot', False, required_type='bool')

    zone = verify_field_variable_definition(
        simulation, zone_vardef, 'forcing zone %r zone definition' % name
    )
    target = verify_field_variable_definition(
        simulation, target_vardef, 'forcing zone %r target' % name
    )

    verify_key(
        'forcing zone type',
        ztype,
        ('MomentumForcing', 'ScalarForcing'),
        'penalty forcing zone %s' % name,
    )
    if ztype == 'MomentumForcing':
        varname = inp.get_value('variable', 'u', required_type='string')
    else:
        varname = inp.get_value('variable', required_type='string')

    fzones.setdefault(varname, []).append(ForcingZone(name, zone, target, penalty))

    if plot:
        # Save zone blending function to a plot file for easy verification
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        pfile = prefix + '_blending_zone_%s.pvd' % name
        zone.rename('beta', 'beta')
        dolfin.File(pfile) << zone


class ForcingZone:
    def __init__(self, name, zone_func, target_func, penalty):
        self.name = name
        self.beta = zone_func
        self.target = target_func
        self.penalty = dolfin.Constant(penalty)
