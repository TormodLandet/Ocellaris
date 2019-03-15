# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import verify_key
from . import register_known_field, KnownField


@register_known_field('MaxField')
class MaxField(KnownField):
    description = 'The maximum of two fields'

    def __init__(self, simulation, field_inp):
        """
        A max field, represents the maximum of two fields
        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a maximum field %r' % self.name)

    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        field_name0 = field_inp.get_value('field0', required_type='string')
        field_name1 = field_inp.get_value('field1', required_type='string')

        verify_key(
            'field0', field_name0, self.simulation.fields, 'MaxField %s' % self.name
        )
        verify_key(
            'field1', field_name1, self.simulation.fields, 'MaxField %s' % self.name
        )
        self.field0 = self.simulation.fields[field_name0]
        self.field1 = self.simulation.fields[field_name1]

    def get_variable(self, name):
        f0 = self.field0.get_variable(name)
        f1 = self.field1.get_variable(name)
        return dolfin.conditional(f0 >= f1, f0, f1)


@register_known_field('MinField')
class MinField(KnownField):
    description = 'The minimum of two fields'

    def __init__(self, simulation, field_inp):
        """
        A max field, represents the minimum of two fields
        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a minimum field %r' % self.name)

    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        field_name0 = field_inp.get_value('field0', required_type='string')
        field_name1 = field_inp.get_value('field1', required_type='string')

        verify_key(
            'field0', field_name0, self.simulation.fields, 'MinField %s' % self.name
        )
        verify_key(
            'field1', field_name1, self.simulation.fields, 'MinField %s' % self.name
        )
        self.field0 = self.simulation.fields[field_name0]
        self.field1 = self.simulation.fields[field_name1]

    def get_variable(self, name):
        f0 = self.field0.get_variable(name)
        f1 = self.field1.get_variable(name)
        return dolfin.conditional(f0 <= f1, f0, f1)
