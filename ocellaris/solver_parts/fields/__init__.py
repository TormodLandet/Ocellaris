# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from ocellaris.utils import ocellaris_error


DEFAULT_POLYDEG = 2
_KNOWN_FIELDS = {}


def add_known_field(name, known_field_class):
    """
    Register a known field
    """
    _KNOWN_FIELDS[name] = known_field_class


def register_known_field(name):
    """
    A class decorator to register known fields
    """

    def register(known_field_class):
        add_known_field(name, known_field_class)
        return known_field_class

    return register


def get_known_field(name):
    """
    Return a known field by name
    """
    try:
        return _KNOWN_FIELDS[name]
    except KeyError:
        ocellaris_error(
            'Field type "%s" not found' % name,
            'Available field types:\n'
            + '\n'.join(
                '  %-20s - %s' % (n, s.description) for n, s in sorted(_KNOWN_FIELDS.items())
            ),
        )
        raise


class KnownField(object):
    description = 'No description available'


from . import scalar_field
from . import vector_field
from . import sharp_field
from . import airy_waves
from . import raschii_waves
from . import wave_outflow
from . import blended_field
from . import conditional_field
from . import free_surface_zone
from . import water_block_field
