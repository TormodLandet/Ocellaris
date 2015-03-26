from . import register_multi_phase_model, MultiPhaseModel
import dolfin

@register_multi_phase_model('SinglePhase')
class SinglePhaseScheme(MultiPhaseModel):
    description = 'A single phase model with a single denisity and laminar viscosity everywhere'
    
    def __init__(self, simulation):
        self.simulation = simulation
        self.rho0 = self.simulation.input.get_value('physical_properties/rho0', 1.0, 'float')
        self.nu0 = self.simulation.input.get_value('physical_properties/nu0', required_type='float')
    
    def get_density(self):
        return dolfin.Constant(self.rho0)

    def get_laminar_kinematic_viscosity(self):
        return dolfin.Constant(self.nu0)
    
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
        mu = self.nu0*self.rho0
        return mu, mu
