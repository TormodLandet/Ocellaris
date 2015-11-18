# encoding: utf-8
from __future__ import division
import dolfin
import numpy
from dolfin import Function, Constant, dot, div, grad, jump, dx, dS
from . import register_multi_phase_model, MultiPhaseModel 
from ..convection import get_convection_scheme


CONVECTION_SCHEME = 'HRIC'
CONTINUOUS_FIELDS = True
CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION = True


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

        # The projected density and viscosity functions for the new time step can be made continuous
        self.continuous_fields = simulation.input.get_value('multiphase_solver/continuous_fields',
                                                            CONTINUOUS_FIELDS, 'bool')
        if self.continuous_fields:
            mesh = simulation.data['mesh']
            P = V.ufl_element().degree()
            V_star = dolfin.FunctionSpace(mesh, 'CG', P+1)
            self.continuous_c_star = dolfin.Function(V_star)
            self.continuous_c = dolfin.Function(V_star)
            self.continuous_c_old = dolfin.Function(V_star)
            
        # Calculate mu from rho and nu (i.e mu is quadratic in c) or directly from c (linear in c)
        self.calculate_mu_directly_from_colour_function = \
            simulation.input.get_value('multiphase_solver/calculate_mu_directly_from_colour_function',
                                       CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION, 'bool')
        
        # Get the physical properties
        self.rho0 = self.simulation.input.get_value('physical_properties/rho0', required_type='float')
        self.rho1 = self.simulation.input.get_value('physical_properties/rho1', required_type='float')
        self.nu0 = self.simulation.input.get_value('physical_properties/nu0', required_type='float')
        self.nu1 = self.simulation.input.get_value('physical_properties/nu1', required_type='float')
        
        # The convection blending function that counteracts numerical diffusion
        scheme = simulation.input.get_value('convection/c/convection_scheme', CONVECTION_SCHEME, 'string')
        scheme_class = get_convection_scheme(scheme)
        self.convection_scheme = scheme_class(simulation, 'c')
        self.convection_scheme_star = scheme_class(simulation, 'c_star')
        self.need_gradient = scheme_class.need_alpha_gradient
        
        # Create the equations when the simulation starts
        self.simulation.hooks.add_pre_simulation_hook(self.on_simulation_start, 'BlendedAlgebraicVofModel setup equations')
        
        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(self.update, 'BlendedAlgebraicVofModel - update colour field')
        
        # Report divergence of the velocity field after each time step
        simulation.hooks.add_post_timestep_hook(self.report_divergence, 'BlendedAlgebraicVofModel - report velocity divergence')

        simulation.log.info('Creating blended VOF multiphase model')
        simulation.log.info('    Using convection scheme %s for the colour function' % scheme)
        simulation.log.info('    Using continuous rho and nu fields: %r' % self.continuous_fields)
    
    def on_simulation_start(self):
        """
        This runs when the simulation starts. It does not run in __init__
        since the solver needs the density and viscosity we define, and
        we need the velocity that is defined by the solver
        """
        beta = self.convection_scheme.blending_function
        
        # The time step (real value to be supplied later)
        self.dt = Constant(1.0)
        
        # Use first order backward time difference on the first time step
        # Coefficients for u, up and upp 
        self.time_coeffs = Constant([1, -1, 0])
        self.extrapolation_coeffs = [1, 0]
        
        # Setup the equation to solve
        c = self.simulation.data['c']
        cp = self.simulation.data['cp']
        cpp = self.simulation.data['cpp']
        c_star = self.simulation.data['c_star']
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('c', [])
        
        # Make sure the convection scheme has something usefull in the first iteration
        c_star.assign(c)
        
        # Define equation for advection of the colour function
        #    ∂c/∂t +  ∇⋅(c u) = 0
        
        # At the previous time step (known advecting velocity "up")   
        vel = self.simulation.data['up']
        self.eq = AdvectionEquation(self.simulation, cp, cpp, vel, beta, self.time_coeffs, dirichlet_bcs)
        
        # At the next time step (extrapolated advecting velocity)
        Vu = self.simulation.data['Vu']
        self.u_conv_comps = [dolfin.Function(Vu) for _ in range(self.simulation.ndim)]
        self.u_conv = dolfin.as_vector(self.u_conv_comps)
        beta_star = self.convection_scheme_star.blending_function
        self.eq_star = AdvectionEquation(self.simulation, c, cp, self.u_conv, beta_star, self.time_coeffs, dirichlet_bcs)
        
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
            if self.continuous_fields:
                c = self.continuous_c
            else:
                c = self.simulation.data['cp']
        elif k == -1:
            if self.continuous_fields:
                c = self.continuous_c_old
            else:
                c = self.simulation.data['cpp']
        elif k == 1:
            if self.continuous_fields:
                c = self.continuous_c_star
            else:
                c = self.simulation.data['c_star']
        
        if force_steady:
            c = self.simulation.data['c']
        
        if force_boundedness:
            c = dolfin.max_value(dolfin.min_value(c, Constant(1.0)), Constant(0.0))
        
        if sharp_interface:
            c = dolfin.conditional(dolfin.ge(c, 0.5), Constant(1.0), Constant(0.0))
        
        return c
    
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
        return Constant(self.rho0)*c + Constant(self.rho1)*(1 - c)
    
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
        c = self.get_colour_function(k)
        return Constant(self.nu0)*c + Constant(self.nu1)*(1 - c)
    
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
            c = self.get_colour_function(k)
            mu0 = self.nu0*self.rho0
            mu1 = self.nu1*self.rho1
            return Constant(mu0)*c + Constant(mu1)*(1 - c)
        
        else:
            nu = self.get_laminar_kinematic_viscosity(k, c)
            rho = self.get_density(k, c)
            return nu*rho
    
    def get_density_range(self):
        """
        Return the maximum and minimum densities, rho
        """
        return min(self.rho0, self.rho1), max(self.rho0, self.rho1) 
               
    def get_laminar_kinematic_viscosity_range(self):
        """
        Return the maximum and minimum kinematic viscosities, nu
        """
        return min(self.nu0, self.nu1), max(self.nu0, self.nu1)
    
    def get_laminar_dynamic_viscosity_range(self):
        """
        The minimum and maximum laminar dynamic viscosities, mu.
        
        Mu is either calculated directly from the colour function, in this
        case mu is a linear function, or as a product of nu and rho, where
        it is a quadratic function and can have (in i.e the case of water
        and air) have maximum value in the middle of the range c ∈ (0, 1)
        """
        if self.calculate_mu_directly_from_colour_function:
            mu0 = self.nu0*self.rho0
            mu1 = self.nu1*self.rho1
            return min(mu0, mu1), max(mu0, mu1)
        else:
            c = numpy.linspace(0, 1, 1000)
            nu = self.nu0*c + self.nu1*(1 - c)
            rho = self.rho0*c + self.rho1*(1 - c)
            mu = nu*rho
            return mu.min(), mu.max()
    
    def update(self, it, t, dt):
        """
        Update the VOF field by advecting it for a time dt
        using the given divergence free velocity field
        """
        timer = dolfin.Timer('Ocellaris update VOF')
        self.dt.assign(dt)
        
        # Reconstruct the gradients
        if self.need_gradient:
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
            A = self.eq.assemble_lhs()
            b = self.eq.assemble_rhs()
            dolfin.solve(A, c.vector(), b)
        
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
        A_star = self.eq_star.assemble_lhs()
        b_star = self.eq_star.assemble_rhs()
        dolfin.solve(A_star, c_star.vector(), b_star)
        
        # Optionally use a continuous predicted colour field
        if self.continuous_fields:
            #self.continuous_c_star.interpolate(c_star)
            #self.continuous_c.interpolate(c)
            #self.continuous_c_old.interpolate(cp)
            Vcg = self.continuous_c_star.function_space()          
            self.continuous_c_star.assign(dolfin.project(c_star, Vcg))
            self.continuous_c.assign(dolfin.project(c, Vcg))
            self.continuous_c_old.assign(dolfin.project(cp, Vcg))
        
        # Report total mass balance and divergence
        sum_c = dolfin.assemble(c*dolfin.dx)
        arr_c = c.vector().get_local()
        min_c = dolfin.MPI.min(dolfin.mpi_comm_world(), float(arr_c.min()))
        max_c = dolfin.MPI.max(dolfin.mpi_comm_world(), float(arr_c.max()))
        self.simulation.reporting.report_timestep_value('sum(c)', sum_c)
        self.simulation.reporting.report_timestep_value('min(c)', min_c)
        self.simulation.reporting.report_timestep_value('max(c)', max_c)
        
        # Update the previous values for the next time step
        cpp.assign(cp)
        cp.assign(c)
        
        # Use second order backward time difference after the first time step
        self.time_coeffs.assign(Constant([3/2, -2, 1/2]))
        self.extrapolation_coeffs = [2, -1]
        
        timer.stop()
    
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
            div_u.vector().abs()
            maxdiv = div_u.vector().max()
            self.simulation.reporting.report_timestep_value('max(div(u)|%s)' % fspace_name, maxdiv)


class AdvectionEquation(object):
    def __init__(self, simulation, cp, cpp, u_conv, beta, time_coeffs, dirichlet_bcs):
        """
        This class assembles the advection equation for the colour function, both CG and DG 
        """
        self.simulation = simulation
        self.cp = cp
        self.cpp = cpp
        self.u_conv = u_conv
        self.beta = beta
        self.time_coeffs = time_coeffs
        self.dirichlet_bcs = dirichlet_bcs
        
        # Discontinuous or continuous elements
        Vc_family = simulation.data['Vc'].ufl_element().family()
        self.colour_is_discontinuous = (Vc_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_advection_equation()
    
    def define_advection_equation(self):
        """
        Setup the advection equation for the colour function
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Trial and test functions
        Vc = sim.data['Vc']
        c = dolfin.TrialFunction(Vc)
        d = dolfin.TestFunction(Vc)
        
        c1, c2, c3 = self.time_coeffs
        dt = sim.data['dt']
        u_conv = self.u_conv
        
        if not self.colour_is_discontinuous:
            # Continous Galerkin implementation of the advection equation
            eq = (c1*c + c2*self.cp + c3*self.cpp)/dt*d*dx + div(c*u_conv)*d*dx
        
        else:
            # Upstream and downstream normal velocities
            flux_nU = c*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
            flux_nD = c*(dot(u_conv, n) - abs(dot(u_conv, n)))/2
            
            # Define the blended flux
            # The blending factor beta is not DG, so beta('+') == beta('-')
            b = self.beta('+')
            flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
            
            # Discontinuous Galerkin implementation of the advection equation 
            eq = (c1*c + c2*self.cp + c3*self.cpp)/dt*d*dx \
                 - dot(c*u_conv, grad(d))*dx \
                 + flux*jump(d)*dS
            
            # Enforce Dirichlet BCs weakly
            for dbc in self.dirichlet_bcs:
                eq += dot(u_conv, n)*d*(c - dbc.func())*dbc.ds()
        
        a, L = dolfin.system(eq)
        self.form_lhs = a
        self.form_rhs = L
        self.tensor_lhs = None
        self.tensor_rhs = None
    
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

