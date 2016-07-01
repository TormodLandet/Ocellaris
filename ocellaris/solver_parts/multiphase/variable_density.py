# encoding: utf-8
from __future__ import division
import dolfin
from dolfin import Function, Constant
from . import register_multi_phase_model, MultiPhaseModel
from ocellaris.solver_parts import SlopeLimiter
from ocellaris.utils import linear_solver_from_input, ocellaris_interpolate
from .advection_equation import AdvectionEquation


# Default values, can be changed in the input file
SOLVER = 'gmres'
PRECONDITIONER = 'default'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'relative_tolerance': 1e-15,
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
        self.use_analytic_solution = simulation.input.get_value('multiphase_solver/analytic_solution', False, 'bool') 
    
    def on_simulation_start(self):
        """
        This runs when the simulation starts. It does not run in __init__
        since the N-S solver needs the density and viscosity we define, and
        we need the velocity that is defined by the solver
        """
        # Use backward euler (BDF1) for timestep 1 
        self.time_coeffs = Constant([1, -1, 0])
        
        if self.use_analytic_solution:
            return
        
        if dolfin.norm(self.rho_pp.vector()) > 0:
            # Use BDF2 from the start
            self.time_coeffs = Constant([3/2, -2, 1/2])
            self.simulation.log.info('Using second order timestepping from the start in VariableDensity')
        
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
    
    def get_density(self, k):
        """
        Return the function as defined on timestep t^{n+k}
        """
        if k == 0:
            return self.rho
        elif k == -1:
            return self.rho_p
        elif k == -2:
            return self.rho_pp
    
    def get_laminar_kinematic_viscosity(self, k):
        """
        It is assumed that the kinematic viscosity is constant
        """
        return Constant(self.nu)
    
    def get_laminar_dynamic_viscosity(self, k):
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
        
        if it != 1:
            # Update the previous values
            self.rho_pp.assign(self.rho_p)
            self.rho_p.assign(self.rho)
            self.time_coeffs.assign(Constant([3/2, -2, 1/2]))
        
        if self.use_analytic_solution:
            # Use the prescribed field for testing
            cpp_code = self.simulation.input.get_value('initial_conditions/rho_p/cpp_code',
                                                       required_type='string')
            description = 'initial condition for rho_p'
            V = self.simulation.data['Vrho']
            ocellaris_interpolate(self.simulation, cpp_code, description, V, self.rho)
        else:
            # Solve the advection equation for the density rho
            A = self.eq.assemble_lhs()
            b = self.eq.assemble_rhs()
            dolfin.solve(A, self.rho.vector(), b)
            self.slope_limiter.run()
        
        self.simulation.reporting.report_timestep_value('min(rho)', self.rho.vector().min())
        self.simulation.reporting.report_timestep_value('max(rho)', self.rho.vector().max())
        
        timer.stop()
