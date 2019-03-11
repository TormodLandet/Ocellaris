# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from . import register_multi_phase_model, MultiPhaseModel
import dolfin


@register_multi_phase_model('SinglePhase')
class SinglePhaseScheme(MultiPhaseModel):
    description = 'A single phase model with a single denisity and laminar viscosity everywhere'

    def __init__(self, simulation):
        self.simulation = simulation
        simulation.log.info('Creating single phase model')
        props = self.simulation.input.get_value('physical_properties', {}, 'Input')
        self.rho0 = props.get_value('rho', 1.0, 'float')
        self.nu0 = props.get_value('nu', None, 'float')
        if self.nu0 is None:
            self.nu0 = 1e100
            simulation.log.warning(
                'No viscosity (physical_properties/nu) given. Using %r' % self.nu0
                + ', make sure to change this if you are using the viscosity!'
            )
        simulation.log.info('    rho = %r' % self.rho0)
        simulation.log.info('    nu = %r' % self.nu0)

    def get_density(self, k):
        """
        Get fluid density function rho at timestep t^{n+k}
        """
        return dolfin.Constant(self.rho0)

    def get_laminar_kinematic_viscosity(self, k):
        """
        Get fluid kinematic viscosity function nu at timestep t^{n+k}
        """
        return dolfin.Constant(self.nu0)

    def get_laminar_dynamic_viscosity(self, k):
        """
        Get fluid dynamic viscosity function mu at timestep t^{n+k}
        """
        return dolfin.Constant(self.nu0 * self.rho0)

    def get_density_range(self):
        """
        The minimum and maximum fluid densities
        These will be identical for single phase flows
        """
        return self.rho0, self.rho0

    def get_laminar_kinematic_viscosity_range(self):
        """
        The minimum and maximum laminar kinematic viscosities
        These will be identical for single phase flows
        """
        return self.nu0, self.nu0

    def get_laminar_dynamic_viscosity_range(self):
        """
        The minimum and maximum laminar dynamic viscosities
        These will be identical for single phase flows
        """
        mu = self.nu0 * self.rho0
        return mu, mu
