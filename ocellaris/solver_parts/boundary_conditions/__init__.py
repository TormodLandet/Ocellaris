# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import ocellaris_error, RunnablePythonString


_BOUNDARY_CONDITIONS = {}


def add_boundary_condition(name, boundary_condition_class):
    """
    Register a boundary condition
    """
    _BOUNDARY_CONDITIONS[name] = boundary_condition_class


def register_boundary_condition(name):
    """
    A class decorator to register boundary conditions
    """

    def register(boundary_condition_class):
        add_boundary_condition(name, boundary_condition_class)
        return boundary_condition_class

    return register


def get_boundary_condition(name):
    """
    Return a boundary condition by name
    """
    try:
        return _BOUNDARY_CONDITIONS[name]
    except KeyError:
        ocellaris_error(
            'Boundary condition "%s" not found' % name,
            'Available boundary conditions:\n'
            + '\n'.join(
                '  %-20s - %s' % (n, s.description)
                for n, s in sorted(_BOUNDARY_CONDITIONS.items())
            ),
        )
        raise


class BoundaryConditionCreator(object):
    description = 'No description available'


from .boundary_region import BoundaryRegion
from .dof_marker import get_dof_region_marks, mark_cell_layers

from . import dirichlet
from . import neumann
from . import robin
from . import slip_length
from . import wall
from . import outlet
from . import decomposed
