# encoding: utf-8
from __future__ import division
import dolfin
from dolfin import Function, Constant, FacetNormal, solve
from . import register_multi_phase_model, MultiPhaseModel 
from ..convection import get_convection_scheme
from ..solvers.ipcs_equations import define_advection_problem

@register_multi_phase_model('BlendedAlgebraicVOF')
class BlendedAlgebraicVofModel(MultiPhaseModel):
    description = 'A blended algebraic VOF scheme implementing HRIC/CICSAM type schemes'
    
    def __init__(self, simulation):
        """
        A blended algebraic VOF scheme works by using a specific 
        convection scheme in the advection of the colour function
        that ensures a sharp interface.
        
        * The convection scheme should be the name of a convection
          scheme that is tailored for advection of the colour 
          function, i.e "HRIC", "MHRIC", "RHRIC" etc, 
        * The velocity field should be divergence free
        
        The colour function is unity when rho=rho0 and nu=nu0 and
        zero when rho=rho1 and nu=nu1
        """
        self.simulation = simulation
        
        # Define function space and solution function
        V = simulation.data['Vc']
        simulation.data['c'] = Function(V)
        simulation.data['cp'] = Function(V)
        simulation.data['cpp'] = Function(V)
        simulation.data['c_star']= Function(V)
        
        # Get the physical properties
        self.rho0 = self.simulation.input.get_value('physical_properties/rho0', required_type='float')
        self.rho1 = self.simulation.input.get_value('physical_properties/rho1', required_type='float')
        self.nu0 = self.simulation.input.get_value('physical_properties/nu0', required_type='float')
        self.nu1 = self.simulation.input.get_value('physical_properties/nu1', required_type='float')
        
        # The convection blending function that counteracts numerical diffusion
        scheme = simulation.input.get_value('convection/c/convection_scheme', 'HRIC', 'string')
        scheme_class = get_convection_scheme(scheme)
        self.convection_scheme = scheme_class(simulation, 'c')
        self.convection_scheme_star = scheme_class(simulation, 'c_star')
        
        self.need_gradient = True
        if scheme == 'Upwind':
            self.need_gradient = False
        
        # Create the equations when the simulation starts
        self.simulation.hooks.add_pre_simulation_hook(self.on_simulation_start, 'BlendedAlgebraicVofModel setup equations')
        
        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(self.update, 'BlendedAlgebraicVofModel - update colour field')
        
        # Report divergence of the velocity field after each time step
        simulation.hooks.add_post_timestep_hook(self.report_divergence, 'BlendedAlgebraicVofModel - report velocity divergence')
    
    def on_simulation_start(self):
        """
        This runs when the simulation starts. It does not run in __init__
        since the solver needs the density and viscosity we define, and
        we need the velocity that is defined by the solver
        """
        mesh = self.simulation.data['mesh']
        beta = self.convection_scheme.blending_function
        
        # The time step (real value to be supplied later)
        self.dt = Constant(1.0)
        
        # Use first order backward time difference on the first time step
        # Coefficients for u, up and upp 
        self.time_coeffs = Constant([1, -1, 0])
        self.extrapolation_coeffs = [1, 0]
        
        # The normal on each face
        normal = FacetNormal(mesh)
        
        # Setup the equation to solve
        V = self.simulation.data['Vc']
        c = self.simulation.data['c']
        cp = self.simulation.data['cp']
        cpp = self.simulation.data['cpp']
        c_star = self.simulation.data['c_star']
        trial = dolfin.TrialFunction(V)
        test = dolfin.TestFunction(V)
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('c', [])
        
        # Make sure the convection scheme has something usefull in the first iteration
        c_star.assign(c)
        
        # Define equation for advection of the colour function
        #    ∂c/∂t +  ∇⋅(c u) = 0
        
        # At the previous time step (known advecting velocity "up")   
        vel = self.simulation.data['up']
        r = Constant(1.0)
        thetas = Constant([1.0, 0.0, 0.0])
        self.eq = define_advection_problem(trial, test, cp, cpp, vel, r, normal, beta,
                                           self.time_coeffs, thetas, self.dt, dirichlet_bcs)
        
        
        # At the next time step (extrapolated advecting velocity)
        Vu = self.simulation.data['Vu']
        self.u_conv_comps = [dolfin.Function(Vu) for _ in range(self.simulation.ndim)]
        self.u_conv = dolfin.as_vector(self.u_conv_comps)
        beta_star = self.convection_scheme_star.blending_function
        self.eq_star = define_advection_problem(trial, test, c, cp, self.u_conv, r, normal, beta_star,
                                                self.time_coeffs, thetas, self.dt, dirichlet_bcs)
        
        # Add some debugging plots to show results in 2D
        self.simulation.plotting.add_plot('c', c, clim=(0, 1))        
        self.simulation.plotting.add_plot('c_beta', beta)
        
        if self.need_gradient:
            # Reconstruct the gradient from the colour function DG0 field
            self.convection_scheme.gradient_reconstructor.initialize()
            gradient = self.convection_scheme.gradient_reconstructor.gradient
            self.simulation.plotting.add_plot('c_grad', gradient)
    
    def get_colour_function(self, k, force_boundedness=False, sharp_interface=False, force_steady=False):
        """
        Return the colour function on timestep t^{n+k}
        """
        if k == 0:
            c = self.simulation.data['cp']
        elif k == -1:
            c = self.simulation.data['cpp']
        elif k == 1:
            c = self.simulation.data['c_star']
        
        if force_steady:
            c = self.simulation.data['c']
        
        if force_boundedness:
            c = dolfin.max_value(dolfin.min_value(c, Constant(1.0)), Constant(0.0))
        
        if sharp_interface:
            c = dolfin.conditional(dolfin.ge(c, 0.5), Constant(1.0), Constant(0.0))
        
        return c
    
    def get_density(self, k):
        """
        Calculate the blended density function as a weighted sum of
        rho0 and rho1. The colour function is unity when rho=rho0
        and zero when rho=rho1
        
        Return the function as defined on timestep t^{n+k}
        """
        c = self.get_colour_function(k)
        return Constant(self.rho0)*c + Constant(self.rho1)*(1 - c)
    
    def get_laminar_kinematic_viscosity(self, k):
        """
        Calculate the blended kinematic viscosity function as a weighted
        sum of nu0 and nu1. The colour function is unity when nu=nu0 and
        zero when nu=nu1
        
        Return the function as defined on timestep t^{n+k}
        """
        c = self.get_colour_function(k)
        return Constant(self.nu0)*c + Constant(self.nu1)*(1 - c)
    
    def get_density_range(self):
        """
        Return the maximum and minimum densities
        """
        return min(self.rho0, self.rho1), max(self.rho0, self.rho1) 
               
    def get_laminar_kinematic_viscosity_range(self):
        """
        Return the maximum and minimum kinematic viscosities
        """
        return min(self.nu0, self.nu1), max(self.nu0, self.nu1)
    
    def get_laminar_dynamic_viscosity_range(self):
        """
        The minimum and maximum laminar dynamic viscosities
        """
        mu0 = self.nu0*self.rho0
        mu1 = self.nu1*self.rho1
        return min(mu0, mu1), max(mu0, mu1)
    
    def update(self, it, t, dt):
        """
        Update the VOF field by advecting it for a time dt
        using the given divergence free velocity field
        """
        self.dt.assign(dt)
        
        if self.need_gradient:
            # Reconstruct the gradients
            self.convection_scheme.gradient_reconstructor.reconstruct()
            self.convection_scheme_star.gradient_reconstructor.reconstruct()
        
        # Update the convection blending factors
        vel = self.simulation.data['up']
        self.convection_scheme.update(t, dt, vel)
        
        # Get the functions
        c = self.simulation.data['c']
        cp = self.simulation.data['cp']
        cpp = self.simulation.data['cpp']
        
        # Solve the advection equations for the colour field
        if it == 1:
            c.assign(self.simulation.data['cp'])
        else: 
            a, L = self.eq
            solve(a == L, c)
        
        # Update the extrapolated convecting velocity
        e1, e2 = self.extrapolation_coeffs
        for d, uci in enumerate(self.u_conv_comps):
            uciv = uci.vector() 
            uciv.zero()
            uciv.axpy(e1, self.simulation.data['up%d' % d].vector())
            uciv.axpy(e2, self.simulation.data['upp%d' % d].vector())
        
        # Update the convection blending factors
        self.convection_scheme_star.update(t, dt, self.u_conv)
        
        # Solve the advection equations for the extrapolated colour field
        c_star = self.simulation.data['c_star']
        a_star, L_star = self.eq_star
        solve(a_star == L_star, c_star)
        
        # Report total mass balance and divergence
        sum_c = dolfin.assemble(c*dolfin.dx)
        arr_c = c.vector().array()
        min_c = arr_c.min()
        max_c = arr_c.max()
        self.simulation.reporting.report_timestep_value('sum(c)', sum_c)
        self.simulation.reporting.report_timestep_value('min(c)', min_c)
        self.simulation.reporting.report_timestep_value('max(c)', max_c)
        
        # Update the previous values for the next time step
        cpp.assign(cp)
        cp.assign(c)
        
        # Use second order backward time difference after the first time step
        self.time_coeffs.assign(Constant([3/2, -2, 1/2]))
        self.extrapolation_coeffs = [2, -1]
    
    def report_divergence(self):
        """
        It is very important in VOF that the velocity field is divergence free.
        We report the divergence in the function space of the colour function
        so that this can be checked
        """
        vel = self.simulation.data['up']
        
        for fspace_name in ('Vc', 'Vu', 'Vp'):
            V = self.simulation.data[fspace_name]
            div_u = dolfin.project(dolfin.nabla_div(vel), V)
            maxdiv = abs(div_u.vector().array()).max()
            self.simulation.reporting.report_timestep_value('max(div(u)|%s)' % fspace_name, maxdiv)
