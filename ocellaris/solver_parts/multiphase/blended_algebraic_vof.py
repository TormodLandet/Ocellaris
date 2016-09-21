# encoding: utf-8
from __future__ import division
import dolfin
from dolfin import Function, Constant
from . import register_multi_phase_model, MultiPhaseModel
from ..convection import get_convection_scheme, StaticScheme
from .vof import VOFMixin
from .advection_equation import AdvectionEquation


# Default values, can be changed in the input file
CONVECTION_SCHEME = 'Upwind'
CONTINUOUS_FIELDS = False
CALCULATE_MU_DIRECTLY_FROM_COLOUR_FUNCTION = False
FORCE_STEADY = False
FORCE_BOUNDED = False
FORCE_SHARP = False


@register_multi_phase_model('BlendedAlgebraicVOF')
class BlendedAlgebraicVofModel(VOFMixin, MultiPhaseModel):
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
        simulation.data['c_star'] = Function(V)

        # The projected density and viscosity functions for the new time step can be made continuous
        self.continuous_fields = simulation.input.get_value('multiphase_solver/continuous_fields',
                                                            CONTINUOUS_FIELDS, 'bool')
        if self.continuous_fields:
            mesh = simulation.data['mesh']
            P = V.ufl_element().degree()
            V_star = dolfin.FunctionSpace(mesh, 'CG', P + 1)
            self.continuous_c_star = dolfin.Function(V_star)
            self.continuous_c = dolfin.Function(V_star)
            self.continuous_c_old = dolfin.Function(V_star)
            
        self.force_steady = simulation.input.get_value('multiphase_solver/force_steady', FORCE_STEADY, 'bool')
        self.force_bounded = simulation.input.get_value('multiphase_solver/force_bounded', FORCE_BOUNDED, 'bool')
        self.force_sharp = simulation.input.get_value('multiphase_solver/force_sharp', FORCE_SHARP, 'bool')
            
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
        simulation.hooks.add_pre_simulation_hook(self.on_simulation_start, 'BlendedAlgebraicVofModel setup equations')
        
        # Update the rho and nu fields before each time step
        simulation.hooks.add_pre_timestep_hook(self.update, 'BlendedAlgebraicVofModel - update colour field')
        
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
        c.assign(cp)
        c_star.assign(cp)
        
        # Define equation for advection of the colour function
        #    ∂c/∂t +  ∇⋅(c u) = 0
        
        # At the previous time step (known advecting velocity "up")
        vel = self.simulation.data['up']
        Vc = self.simulation.data['Vc']
        self.eq = AdvectionEquation(self.simulation, Vc, cp, cpp, vel, beta, self.time_coeffs, dirichlet_bcs)
        
        # At the next time step (extrapolated advecting velocity)
        Vu = self.simulation.data['Vu']
        self.u_conv_comps = [dolfin.Function(Vu) for _ in range(self.simulation.ndim)]
        self.u_conv = dolfin.as_vector(self.u_conv_comps)
        beta_star = self.convection_scheme_star.blending_function
        self.eq_star = AdvectionEquation(self.simulation, self.simulation.data['Vc'],
                                         c, cp, self.u_conv, beta_star,
                                         self.time_coeffs, dirichlet_bcs)
        
        # Add some debugging plots to show results in 2D
        self.simulation.plotting.add_plot('c', c, clim=(0, 1))
        self.simulation.plotting.add_plot('c_beta', beta)
        
        if self.need_gradient:
            # Reconstruct the gradient from the colour function DG0 field
            self.convection_scheme.gradient_reconstructor.initialize()
            gradient = self.convection_scheme.gradient_reconstructor.gradient
            self.simulation.plotting.add_plot('c_grad', gradient)
    
    def get_colour_function(self, k):
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
        
        if self.force_steady:
            c = self.simulation.data['c']
        
        if self.force_bounded:
            c = dolfin.max_value(dolfin.min_value(c, Constant(1.0)), Constant(0.0))
        
        if self.force_sharp:
            c = dolfin.conditional(dolfin.ge(c, 0.5), Constant(1.0), Constant(0.0))
        
        return c
    
    def update(self, timestep_number, t, dt):
        """
        Update the VOF field by advecting it for a time dt
        using the given divergence free velocity field
        """
        timer = dolfin.Timer('Ocellaris update VOF')
        self.dt.assign(dt)
        is_static = isinstance(self.convection_scheme, StaticScheme)
        
        # Reconstruct the gradients
        if self.need_gradient:
            self.convection_scheme.gradient_reconstructor.reconstruct()
            self.convection_scheme_star.gradient_reconstructor.reconstruct()
        
        # Update the convection blending factors
        if not is_static:
            vel = self.simulation.data['up']
            self.convection_scheme.update(t, dt, vel)
        
        # Get the functions
        c = self.simulation.data['c']
        cp = self.simulation.data['cp']
        cpp = self.simulation.data['cpp']
        c_star = self.simulation.data['c_star']
        
        # Solve the advection equations for the colour field
        if timestep_number == 1 or is_static:
            c.assign(cp)
        else:
            A = self.eq.assemble_lhs()
            b = self.eq.assemble_rhs()
            dolfin.solve(A, c.vector(), b)
        
        if is_static:
            c_star.assign(cp)
        else:
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
            A_star = self.eq_star.assemble_lhs()
            b_star = self.eq_star.assemble_rhs()
            dolfin.solve(A_star, c_star.vector(), b_star)
        
        # Optionally use a continuous predicted colour field
        if self.continuous_fields:
            Vcg = self.continuous_c.function_space()
            dolfin.project(c_star, Vcg, function=self.continuous_c_star)
            dolfin.project(c, Vcg, function=self.continuous_c)
            dolfin.project(cp, Vcg, function=self.continuous_c_old)
        
        # Report properties of the colour field
        sum_c = dolfin.assemble(c * dolfin.dx)
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
