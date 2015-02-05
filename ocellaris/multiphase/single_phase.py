from . import register_multi_phase_model, MultiPhaseModel
import dolfin

@register_multi_phase_model('SinglePhase')
class SinglePhaseScheme(MultiPhaseModel):
    description = 'A single phase model with a single denisity and laminar viscosity everywhere'
    
    def __init__(self, simulation):
        self.simulation = simulation
    
    def get_density(self):
        rho0 = self.simulation.input['fluid_properties']['rho0']
        return dolfin.Constant(rho0)

    def get_laminar_kinematic_viscosity(self):
        nu0 = self.simulation.input['fluid_properties']['nu0']
        return dolfin.Constant(nu0)

