# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import ocellaris_error


DEFAULT_FAMILY_U = 'DG'
DEFAULT_FAMILY_P = 'DG'
DEFAULT_DEGREE_U = 2
DEFAULT_DEGREE_P = 1

_SOLVERS = {}


def add_solver(name, solver_class):
    """
    Register a Navier-Stokes solver
    """
    _SOLVERS[name] = solver_class


def register_solver(name):
    """
    A class decorator to register Navier-Stokes solvers
    """

    def register(solver_class):
        add_solver(name, solver_class)
        return solver_class

    return register


def get_solver(name):
    """
    Return a Navier-Stokes solver by name
    """
    try:
        return _SOLVERS[name]
    except KeyError:
        ocellaris_error(
            'Navier-Stokes solver "%s" not found' % name,
            'Available solvers:\n'
            + '\n'.join(
                '  %-20s - %s' % (n, s.description) for n, s in sorted(_SOLVERS.items())
            ),
        )
        raise


class Solver(object):
    description = 'No description available'

    @classmethod
    def create_function_spaces(cls, simulation):
        """
        Function space setup for standard flow solvers
        """
        # Get function space names
        Vu_name = simulation.input.get_value(
            'solver/function_space_velocity', DEFAULT_FAMILY_U, 'string'
        )
        Vp_name = simulation.input.get_value(
            'solver/function_space_pressure', DEFAULT_FAMILY_P, 'string'
        )

        # Get function space polynomial degrees
        Pu = simulation.input.get_value(
            'solver/polynomial_degree_velocity', DEFAULT_DEGREE_U, 'int'
        )
        Pp = simulation.input.get_value(
            'solver/polynomial_degree_pressure', DEFAULT_DEGREE_P, 'int'
        )

        # Get the constrained domain
        cd = simulation.data['constrained_domain']
        if cd is None:
            simulation.log.info(
                'Creating function spaces without periodic boundaries (none found)'
            )
        else:
            simulation.log.info('Creating function spaces with periodic boundaries')

        # Create the Navier-Stokes function spaces
        mesh = simulation.data['mesh']
        Vu = dolfin.FunctionSpace(mesh, Vu_name, Pu, constrained_domain=cd)
        Vp = dolfin.FunctionSpace(mesh, Vp_name, Pp, constrained_domain=cd)
        simulation.data['Vu'] = Vu
        simulation.data['Vp'] = Vp


class BaseEquation(object):
    # Will be shadowed by object properties after first assemble
    tensor_lhs = None
    tensor_rhs = None

    def assemble_lhs(self):
        if self.tensor_lhs is None:
            self.tensor_lhs = dolfin.assemble(self.form_lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            self.tensor_rhs = dolfin.assemble(self.form_rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs


# Timestepping methods
BDF = 'BDF'
CRANK_NICOLSON = 'CN'


# Flux types
BLENDED = 'Blended'
UPWIND = 'Upwind'
LOCAL_LAX_FRIEDRICH = 'Local Lax-Friedrich'


# Velocity post-processing
BDM = 'BDM'


from . import analytical_solution
from . import coupled
from . import coupled_ldg
from . import ipcs
from . import ipcs_algebraic
from . import simple
from . import pimple
