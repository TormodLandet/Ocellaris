# encoding: utf-8
from __future__ import division
import dolfin
from dolfin import Function, Constant
from . import register_multi_phase_model, MultiPhaseModel
from ocellaris.solver_parts import SlopeLimiter
from ocellaris.utils import linear_solver_from_input
from .advection_equation import AdvectionEquation


# Default values, can be changed in the input file
SOLVER = 'gmres'
PRECONDITIONER = 'default'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'relative_tolerance': 1e-10,
                     'absolute_tolerance': 1e-15}


@register_multi_phase_model('VariableDensity')
class VariableDensityModel(MultiPhaseModel):
    description = 'A variable density DG transport equation model'
    
    def __init__(self, simulation):
        """
        A variable density scheme solving D/Dt(rho)=0 with a
        constant kinematic vosicsity nu and a variable dynamic
        visocisty mu (depending on rho)        
        """
        self.simulation = simulation
        simulation.log.info('Creating variable density multiphase model')
        
        # Define function space and solution function
        V = simulation.data['Vrho']
        self.rho = simulation.data['rho'] = Function(V)
        self.rho_p = simulation.data['rho_p'] = Function(V)
        self.rho_pp = simulation.data['rho_pp'] = Function(V)
        
        # Get the physical properties
        self.rho_min = self.simulation.input.get_value('physical_properties/rho_min', required_type='float')
        self.rho_max = self.simulation.input.get_value('physical_properties/rho_max', required_type='float')
        self.nu = self.simulation.input.get_value('physical_properties/nu', required_type='float')
        
        # Create the equations when the simulation starts
        self.simulation.hooks.add_pre_simulation_hook(self.on_simulation_start, 'VariableDensityModel setup equations')
        
        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(self.update, 'VariableDensityModel - update density field')
        
        self.slope_limiter = SlopeLimiter(simulation, 'rho', self.rho)
    
    def on_simulation_start(self):
        """
        This runs when the simulation starts. It does not run in __init__
        since the solver needs the density and viscosity we define, and
        we need the velocity that is defined by the solver
        """
        # The time step (real value to be supplied later)
        self.dt = Constant(1.0)
        
        # Use first order backward time difference on the first time step
        # Coefficients for u, up and upp 
        self.time_coeffs = Constant([1, -1, 0])
        
        # Make sure the convection scheme has something usefull in the first iteration
        self.rho.assign(self.rho_p)
        
        # Define equation for advection of the density
        #    ∂ρ/∂t +  ∇⋅(ρ u) = 0   
        vel = self.simulation.data['up']
        beta = None
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('rho', [])
        self.eq = AdvectionEquation(self.simulation, self.simulation.data['Vrho'],
                                    self.rho_p, self.rho_pp, vel, beta, 
                                    self.time_coeffs, dirichlet_bcs)
        
        self.solver = linear_solver_from_input(self.simulation, 'solver/rho', SOLVER, PRECONDITIONER, None, KRYLOV_PARAMETERS)
        
        # Add some debugging plots to show results in 2D
        self.simulation.plotting.add_plot('rho', self.rho, clim=(self.rho_min, self.rho_max))        
    
    def get_density(self, k=0):
        """
        Return the function as defined on timestep t^{n+k}
        """
        if k == 0:
            return self.rho
        elif k == -1:
            return self.rho_p
        elif k == -2:
            return self.rho_pp
    
    def get_laminar_kinematic_viscosity(self, k=0):
        """
        It is assumed that the kinematic viscosity is constant
        """
        return Constant(self.nu)
    
    def get_laminar_dynamic_viscosity(self, k=0):
        """
        Calculate the blended dynamic viscosity function as a function
        of the (constant) nu and (variable) rho
        
        Return the function as defined on timestep t^{n+k}
        """
        nu = self.get_laminar_kinematic_viscosity(k)
        rho = self.get_density(k)
        return nu*rho
    
    def get_density_range(self):
        """
        Return the maximum and minimum densities, rho
        """
        return self.rho_min, self.rho_max
               
    def get_laminar_kinematic_viscosity_range(self):
        """
        Return the maximum and minimum kinematic viscosities, nu
        """
        return self.nu, self.nu 
    
    def get_laminar_dynamic_viscosity_range(self):
        """
        The minimum and maximum laminar dynamic viscosities, mu.
        """
        return self.nu*self.rho_min, self.nu*self.rho_max
    
    def update(self, it, t, dt):
        """
        Update the density field by advecting it for a time dt
        using the given divergence free velocity field
        """
        timer = dolfin.Timer('Ocellaris update rho')
        self.dt.assign(dt)
        
        if it != 1:
            # Update the previous values
            self.rho_pp.assign(self.rho_p)
            self.rho_p.assign(self.rho)
            self.time_coeffs.assign(Constant([3/2, -2, 1/2]))
        
        A = self.eq.assemble_lhs()
        b = self.eq.assemble_rhs()
        self.solver.solve(A, self.rho.vector(), b)
        self.slope_limiter.run()
        
        self.simulation.reporting.report_timestep_value('min(rho)', self.rho.vector().min())
        self.simulation.reporting.report_timestep_value('max(rho)', self.rho.vector().max())
        
        timer.stop()
