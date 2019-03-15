# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from ocellaris.utils import verify_key, verify_field_variable_definition
from . import register_known_field, KnownField


@register_known_field('BlendedField')
class BlendedField(KnownField):
    description = 'Blended field'

    def __init__(self, simulation, field_inp):
        """
        A blended field

        All variables, such as phi are expressed as

            phi = (1 - b) * phi0 + b * phi1
            b is the scalar blending function with values [0, 1]

        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a blended field %r' % self.name)

    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        field_name0 = field_inp.get_value('field0', required_type='string')
        field_name1 = field_inp.get_value('field1', required_type='string')
        blend_def = field_inp.get_value('blending_function', required_type='string')

        verify_key(
            'blending field0',
            field_name0,
            self.simulation.fields,
            'BlendedField %s' % self.name,
        )
        verify_key(
            'blending field1',
            field_name1,
            self.simulation.fields,
            'BlendedField %s' % self.name,
        )
        blend = verify_field_variable_definition(
            self.simulation, blend_def, 'BlendedField %s' % self.name
        )

        self.field0 = self.simulation.fields[field_name0]
        self.field1 = self.simulation.fields[field_name1]
        self.blending_function = blend

    def get_variable(self, name):
        f0 = self.field0.get_variable(name)
        f1 = self.field1.get_variable(name)
        b = self.blending_function
        return (1 - b) * f0 + b * f1
