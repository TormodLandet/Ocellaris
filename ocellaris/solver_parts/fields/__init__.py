import dolfin
from ocellaris.utils import ocellaris_error


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
        ocellaris_error('Field type "%s" not found' % name,
                        'Available field types:\n' +
                        '\n'.join('  %-20s - %s' % (n, s.description)
                                  for n, s in sorted(_KNOWN_FIELDS.items())))
        raise


class KnownField(object):
    description = 'No description available'


from . import airy_waves