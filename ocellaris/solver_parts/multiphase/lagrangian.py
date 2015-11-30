# encoding: utf-8
from __future__ import division
import dolfin
from . import register_multi_phase_model, MultiPhaseModel
from .vof import VOFMixin


# Default values, can be changed in the input file
CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION = False


@register_multi_phase_model('Lagrangian')
class LagrangianMeshMorpher(VOFMixin, MultiPhaseModel):
    description = 'A purely Lagrangian multiphase model'
    
    def __init__(self, simulation):
        """
        A purely Lagrangian multiphase model. The mesh is moved
        according to the calculated fluid velocity after each time 
        step. This will obviously distort the mesh in allmost all
        calculations.
        
        This was implemented as a stepping stone to ALE, and to test
        hydrostatic pressure calculations where the correct answer
        is zero velocity everywhere for all time and ALE should not
        be necessary. 
        
        To initialise the multi phase field the colour function must
        be specified in the input file (as initial condition for "cp").
        The colour function is unity when rho=rho0 and nu=nu0 and
        zero when rho=rho1 and nu=nu1
        """
        self.simulation = simulation
        
        # Define colour function
        V = simulation.data['Vc']
        c = dolfin.Function(V)
        simulation.data['c'] = c
        simulation.data['cp'] = c # Initial conditions are specified for cp, so we need this alias
        
        # The mesh velocity components
        Vu = simulation.data['Vu']
        u_mesh = []
        for d in range(simulation.ndim):
            umi = dolfin.Function(Vu)
            simulation.data['u_mesh%d' % d] = umi
            u_mesh.append(umi)
        simulation.data['u_mesh'] = dolfin.as_vector(u_mesh)
        
        # Displacement vector function to move the mesh
        mesh = simulation.data['mesh']
        Vd = dolfin.VectorFunctionSpace(mesh, Vu.ufl_element().family(), Vu.ufl_element().degree())
        self.displacement = dolfin.Function(Vd)
        self.assigners = [dolfin.FunctionAssigner(Vd.sub(d), Vu) for d in range(simulation.ndim)]
        
        # Calculate mu from rho and nu (i.e mu is quadratic in c) or directly from c (linear in c)
        self.calculate_mu_directly_from_colour_function = \
            simulation.input.get_value('multiphase_solver/calculate_mu_directly_from_colour_function',
                                       CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION, 'bool')
        
        # Get the physical properties
        self.rho0 = self.simulation.input.get_value('physical_properties/rho0', required_type='float')
        self.rho1 = self.simulation.input.get_value('physical_properties/rho1', required_type='float')
        self.nu0 = self.simulation.input.get_value('physical_properties/nu0', required_type='float')
        self.nu1 = self.simulation.input.get_value('physical_properties/nu1', required_type='float')
        
        # Update the rho and nu fields after each time step
        simulation.hooks.add_post_timestep_hook(self.update, 'LagrangianMeshMorpher - update mesh')
    
        simulation.log.info('Creating Lagrangian mesh morphing multiphase model')
    
    def get_colour_function(self, k):
        """
        The colour function follows the cells and does not ever change
        """
        return self.simulation.data['c']
    
    def update(self):
        """
        Update the mesh position according to the calculated fluid velocities
        """
        timer = dolfin.Timer('Ocellaris move mesh')
        sim = self.simulation
        
        # Get updated mesh velocity (use the fluid velocity)
        for d in range(sim.ndim):
            ui = sim.data['u%d' % d]
            umi = sim.data['u_mesh%d' % d]
            umi.assign(ui)
            self.assigners[d].assign(self.displacement.sub(d), umi)
        self.displacement.vector()[:] *= sim.dt
        
        # Move the mesh according to the mesh velocity
        mesh = sim.data['mesh']
        mesh.move(self.displacement)
        sim.update_mesh_data(connectivity_changed=False)
        
        # Report properties of the colour field
        sum_c = dolfin.assemble(sim.data['c']*dolfin.dx)
        sim.reporting.report_timestep_value('sum(c)', sum_c)
        
        timer.stop()
