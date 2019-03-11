# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
from dolfin import Constant, FunctionSpace
from .level_set_view import LevelSetView


class VOFMixin(object):
    """
    This is a mixin class to avoid having duplicates of the methods calculating
    rho, nu and mu. Any subclass using this mixin must define the method
    "get_colour_function(k)" and can also redefine the boolean property that
    controls the way mu is calculated, "calculate_mu_directly_from_colour_function".
    """

    calculate_mu_directly_from_colour_function = True
    default_polynomial_degree_colour = 0
    _level_set_view = None

    @classmethod
    def create_function_space(cls, simulation):
        mesh = simulation.data['mesh']
        cd = simulation.data['constrained_domain']
        Vc_name = simulation.input.get_value(
            'multiphase_solver/function_space_colour', 'Discontinuous Lagrange', 'string'
        )
        Pc = simulation.input.get_value(
            'multiphase_solver/polynomial_degree_colour',
            cls.default_polynomial_degree_colour,
            'int',
        )
        Vc = FunctionSpace(mesh, Vc_name, Pc, constrained_domain=cd)
        simulation.data['Vc'] = Vc
        simulation.ndofs += Vc.dim()

    def set_physical_properties(self, rho0=None, rho1=None, nu0=None, nu1=None, read_input=False):
        """
        Set rho and nu (density and kinematic viscosity) in both domain 0
        and 1. Either specify all of rho0, rho1, nu0 and nu1 or set
        read_input to True which will read from the physical_properties
        section of the simulation input object.
        """
        sim = self.simulation
        if read_input:
            rho0 = sim.input.get_value('physical_properties/rho0', required_type='float')
            rho1 = sim.input.get_value('physical_properties/rho1', required_type='float')
            nu0 = sim.input.get_value('physical_properties/nu0', required_type='float')
            nu1 = sim.input.get_value('physical_properties/nu1', required_type='float')
        self.df_rho0 = Constant(rho0)
        self.df_rho1 = Constant(rho1)
        self.df_nu0 = Constant(nu0)
        self.df_nu1 = Constant(nu1)
        self.df_smallest_rho = self.df_rho0 if rho0 <= rho1 else self.df_rho1

    def set_rho_min(self, rho_min):
        """
        This is used to bring rho_min closer to rho_max for the initial
        linear solver iterations (to speed up convergence)
        """
        self.df_smallest_rho.assign(Constant(rho_min))

    def get_colour_function(self, k):
        """
        Return the colour function on timestep t^{n+k}
        """
        raise NotImplementedError('The get_colour_function method must be implemented by subclass!')

    def get_density(self, k=None, c=None):
        """
        Calculate the blended density function as a weighted sum of
        rho0 and rho1. The colour function is unity when rho=rho0
        and zero when rho=rho1

        Return the function as defined on timestep t^{n+k}
        """
        if c is None:
            assert k is not None
            c = self.get_colour_function(k)
        else:
            assert k is None
        return self.df_rho0 * c + self.df_rho1 * (1 - c)

    def get_laminar_kinematic_viscosity(self, k=None, c=None):
        """
        Calculate the blended kinematic viscosity function as a weighted
        sum of nu0 and nu1. The colour function is unity when nu=nu0 and
        zero when nu=nu1

        Return the function as defined on timestep t^{n+k}
        """
        if c is None:
            assert k is not None
            c = self.get_colour_function(k)
        else:
            assert k is None
        return self.df_nu0 * c + self.df_nu1 * (1 - c)

    def get_laminar_dynamic_viscosity(self, k=None, c=None):
        """
        Calculate the blended dynamic viscosity function as a weighted
        sum of mu0 and mu1. The colour function is unity when mu=mu0 and
        zero when mu=mu1

        Return the function as defined on timestep t^{n+k}
        """
        if self.calculate_mu_directly_from_colour_function:
            if c is None:
                assert k is not None
                c = self.get_colour_function(k)
            else:
                assert k is None
            mu0 = self.df_nu0 * self.df_rho0
            mu1 = self.df_nu1 * self.df_rho1
            return mu0 * c + mu1 * (1 - c)

        else:
            nu = self.get_laminar_kinematic_viscosity(k, c)
            rho = self.get_density(k, c)
            return nu * rho

    def get_density_range(self):
        """
        Return the maximum and minimum densities, rho
        """
        rho0 = self.df_rho0.values()[0]
        rho1 = self.df_rho1.values()[0]
        return min(rho0, rho1), max(rho0, rho1)

    def get_laminar_kinematic_viscosity_range(self):
        """
        Return the maximum and minimum kinematic viscosities, nu
        """
        nu0 = self.df_nu0.values()[0]
        nu1 = self.df_nu1.values()[0]
        return min(nu0, nu1), max(nu0, nu1)

    def get_laminar_dynamic_viscosity_range(self):
        """
        The minimum and maximum laminar dynamic viscosities, mu.

        Mu is either calculated directly from the colour function, in this
        case mu is a linear function, or as a product of nu and rho, where
        it is a quadratic function and can have (in i.e the case of water
        and air) have maximum values() in the middle of the range c ∈ (0, 1)
        """
        rho0 = self.df_rho0.values()[0]
        rho1 = self.df_rho1.values()[0]
        nu0 = self.df_nu0.values()[0]
        nu1 = self.df_nu1.values()[0]

        if self.calculate_mu_directly_from_colour_function:
            mu0 = nu0 * rho0
            mu1 = nu1 * rho1
            return min(mu0, mu1), max(mu0, mu1)
        else:
            c = numpy.linspace(0, 1, 1000)
            nu = nu0 * c + nu1 * (1 - c)
            rho = rho0 * c + rho1 * (1 - c)
            mu = nu * rho
            return mu.min(), mu.max()

    def get_level_set_view(self):
        """
        Get a view of this VOF field as a level set function
        """
        if self._level_set_view is None:
            self._level_set_view = LevelSetView(self.simulation)
            c = self.get_colour_function(0)
            self._level_set_view.set_density_field(c)
        return self._level_set_view
