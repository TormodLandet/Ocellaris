# encoding: utf-8
from __future__ import division
import numpy
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
        self.mesh = simulation.data['mesh']
        
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
        scheme = get_convection_scheme(simulation.input.get_value('convection/c/convection_scheme', 'HRIC', 'string'))
        self.convection_scheme = scheme(simulation, 'c')
        self.convection_scheme_star = scheme(simulation, 'c_star')
        
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
        beta = self.convection_scheme.blending_function
        
        # The time step (real value to be supplied later)
        self.dt = Constant(1.0)
        
        # Use first order backward time difference on the first time step
        # Coefficients for u, up and upp 
        self.time_coeffs = Constant([1, -1, 0])
        self.extrapolation_coeffs = [1, 0]
        
        # The normal on each face
        normal = FacetNormal(self.mesh)
        
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
        if self.simulation.ndim == 2:
            u_conv0 = dolfin.Function(Vu)
            u_conv1 = dolfin.Function(Vu)
            self.u_conv_comps = [u_conv0, u_conv1]
            self.u_conv = dolfin.as_vector(self.u_conv_comps)
        self.eq_star = define_advection_problem(trial, test, c, cp, self.u_conv, r, normal, beta,
                                                self.time_coeffs, thetas, self.dt, dirichlet_bcs)
        
        # Add some debugging plots to show results in 2D
        self.div_u_Vc = dolfin.Function(V)
        self.simulation.plotting.add_plot('c', c, clim=(0, 1))
        self.simulation.plotting.add_plot('div(u)|Vc', self.div_u_Vc)
        
        self.simulation.plotting.add_plot('c_beta', beta)
        
        if self.need_gradient:
            # Reconstruct the gradient from the colour function DG0 field
            self.convection_scheme.gradient_reconstructor.initialize()
            gradient = self.convection_scheme.gradient_reconstructor.gradient
            self.simulation.plotting.add_plot('c_grad', gradient)
        
        # Move this!
        self.umag = dolfin.Function(V)
        self.simulation.plotting.add_plot('p', self.simulation.data['p'])
        self.simulation.plotting.add_plot('mag(u)', self.umag)
        
    def get_colour_function(self, k, force_boundedness=True, sharp_interface=False, force_steady=False):
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
        if self.need_gradient:
            # Reconstruct the gradients
            self.convection_scheme.gradient_reconstructor.reconstruct()
            self.convection_scheme_star.gradient_reconstructor.reconstruct()
        
        # Update the extrapolated convecting velocity
        e1, e2 = self.extrapolation_coeffs
        for d, uci in enumerate(self.u_conv_comps):
            uciv = uci.vector() 
            uciv.zero()
            uciv.axpy(e1, self.simulation.data['up%d' % d].vector())
            uciv.axpy(e2, self.simulation.data['up%d' % d].vector())
        
        # Update the convection blending factors
        vel = self.simulation.data['up']
        self.convection_scheme.update(t, dt, vel)
        self.convection_scheme_star.update(t, dt, self.u_conv)
        
        # Solve the advection equations for the colour field
        self.dt.assign(dt)
        c = self.simulation.data['c'] 
        a, L = self.eq
        solve(a == L, c)
        
        # Solve the advection equations for the extrapolated colour field
        c_star = self.simulation.data['c_star']
        a, L = self.eq_star
        solve(a == L, c_star)
        
        # Ensure boundedness
        for field in (c_star, c):
            break
            a = field.vector().array()
            a[a>1] = 1
            a[a<0] = 0
            field.vector()[:] = a
        
        # Report total mass balance and divergence
        dolfin.project(dolfin.nabla_div(vel), self.simulation.data['Vc'], function=self.div_u_Vc)
        dolfin.project(dolfin.sqrt(dolfin.dot(vel, vel)), self.simulation.data['Vc'], function=self.umag)
        maxdiv = abs(self.div_u_Vc.vector().array()).max()
        sum_c = dolfin.assemble(c*dolfin.dx)
        arr_c = c.vector().array()
        min_c = arr_c.min()
        max_c = arr_c.max()
        self.simulation.reporting.report_timestep_value('max(div(u_adv_c))', maxdiv)
        self.simulation.reporting.report_timestep_value('sum(c)', sum_c)
        self.simulation.reporting.report_timestep_value('min(c)', min_c)
        self.simulation.reporting.report_timestep_value('max(c)', max_c)
        
        # Update the previous values for the next time step
        cp = self.simulation.data['cp']
        cpp = self.simulation.data['cpp']
        cpp.assign(cp)
        cp.assign(c)
        
        # Use second order backward time difference after the first time step
        self.time_coeffs.assign(Constant([3/2, -2, 1/2]))
        self.extrapolation_coeffs = [2, -1]
        
        # Compress interface
        # Not yet tested with FEniCS 1.5!
        # self.compress(t, dt)
    
    def divergence_free_velocity(self):
        """
        Make the advecting velocity field divergence free in the Vc space
        
        This is unfinished, untested code!
        """
        return self.simulation.data['up'] 
    
        from dolfin import nabla_div, nabla_grad, dot, dx
        
        # Get the advecting velocity field in CR1
        vel = self.advecting_velocity
        dolfin.project(self.simulation.data['up'], self.Vadv, function=vel)
        
        # Define equation that makes the advecting velocity divergence free
        # by using the Helmholtz-Hodge decomposition
        Vc = self.simulation.data['Vc']
        u = dolfin.TrialFunction(self.Vh)
        v = dolfin.TestFunction(self.Vh)
        div_u = dolfin.Function(Vc)
        
        eq = dot(nabla_grad(u), nabla_grad(v))*dx - div_u*v*dx
        a, L = dolfin.system(eq)
        
        phi = Function(self.Vh)
        solver = dolfin.KrylovSolver('minres', 'hypre_amg')
        
        # Remove null space from solution of phi
        null_vec = dolfin.Vector(phi.vector())
        null_vec[:] = 1
        null_vec *= 1/null_vec.norm("l2")
        null_space = dolfin.VectorSpaceBasis([null_vec])
        solver.set_nullspace(null_space)
        
        for i in range(10):
            dolfin.project(nabla_div(vel), Vc, function=div_u)
            
            A, b = dolfin.assemble_system(a, L)
            null_space.orthogonalize(b)
            
            # Update the velocity field
            if self.simulation.timestep == 1 and abs(div_u.vector().array()).max() == 0:
                phi.vector().zero()
            else:
                solver.solve(A, phi.vector(), b)
            
            correction = dolfin.project(nabla_grad(phi), self.Vadv)
            vel.vector().axpy(0.7, correction.vector())
            
            print 'Helmholtz-Hodge ts %5d it %3d divmax %15.3e' % (self.simulation.timestep, i+1, 
                                                                   abs(div_u.vector().array()).max())
        
        return vel
    
    def report_divergence(self):
        """
        It is vitally important that the velocity field is divergence free.
        We report the divergence in the function space of the colour function
        so that this can be checked
        """
        vel = self.simulation.data['u']
        
        for fspace_name in ('Vc', 'Vu', 'Vp'):
            V = self.simulation.data[fspace_name]
            div_u = dolfin.project(dolfin.nabla_div(vel), V)
            maxdiv = abs(div_u.vector().array()).max()
            self.simulation.reporting.report_timestep_value('max(div(u)|%s)' % fspace_name, maxdiv)
    
    def compress(self, t, dt):
        """
        Explicit compression
        """
        #self.simulation.plotting.plot('c', '_uncompr')
        
        compr_fac = self.simulation.input.get_value('multiphase_solver/compression_factor', 0.0, 'float')
        if compr_fac == 0:
            return
        
        facet_info = self.simulation.data['facet_info']
        #cell_info = self.simulation.data['cell_info']
        conFC = self.simulation.data['connectivity_FC']
        ndim = self.simulation.ndim
        
        colour_func_arr = self.colour_function.vector().get_local()
        colour_func_dofmap = self.convection_scheme.alpha_dofmap
        
        gradient_vec = self.convection_scheme.gradient_reconstructor.gradient.vector()
        gradient_dofmap0 = self.convection_scheme.gradient_reconstructor.gradient_dofmap0
        gradient_dofmap1 = self.convection_scheme.gradient_reconstructor.gradient_dofmap1
        
        # Fast functions for length of vectors
        if ndim == 2:
            norm = lambda vec: (vec[0]**2 + vec[1]**2)**0.5
        else:
            norm = lambda vec: (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
            
        # Get a numpy array from a location in a VectorFunctionSpace of the gradient
        gdofs  = (gradient_dofmap0, gradient_dofmap1)
        gradient2array = lambda vecfun, i: numpy.array([vecfun[dm[i]] for dm in gdofs], float)
        
        # Find faces where there is a change of colour between the cells
        EPS = 1e-6
        faces_to_flux = []
        for facet in dolfin.facets(self.mesh):
            fidx = facet.index()
            finfo = facet_info[fidx]
            
            # Find the local cells (the two cells sharing this face)
            connected_cells = conFC(fidx)
            
            if len(connected_cells) != 2:
                # Skip boundary facets
                continue
            
            # Indices of the two local cells
            i0, i1 = connected_cells
            
            # Find colour function in cells 0 and 1
            c0 = colour_func_arr[colour_func_dofmap[i0]]
            c1 = colour_func_arr[colour_func_dofmap[i1]]
            
            if abs(c0 - c1) < EPS:
                # Skip areas of constant colour
                continue
            
            # Facet midpoint
            #face_mp = facet_info[fidx].midpoint
            
            # Velocity at the midpoint (do not care which side of the face)
            #ump = numpy.zeros(2, float)
            #self.velocity_field.eval(ump, face_mp)
            
            # Find a normal pointing out from local cell 0
            normal = finfo.normal
            
            # Find average gradient  
            gradient = 0.5*(gradient2array(gradient_vec, i0) + gradient2array(gradient_vec, i1))
            
            # Find indices of downstream ("D") cell and central ("C") cell
            if numpy.dot(normal, gradient) > 0:
                iC, iD = i0, i1
                cC, cD = c0, c1
            else:
                iC, iD = i1, i0
                cC, cD = c1, c0
            
            # We must allow some "diffusion", otherwise the front will not move
            if cD > 0.9:
                continue
            
            # The colour function values are for sorting purposes only,
            # the others are to calculate the compressive flux on this face
            faces_to_flux.append((cD, cC, iD, iC, gradient, normal, finfo.area))
        
        # Sort to bring the largest values of the colour function in
        # the recipient cell first
        faces_to_flux.sort(reverse=False)
        
        for _, _, iD, iC, gradient, normal, area in faces_to_flux:
            # Find updated colour function in D and C cells
            cC = colour_func_arr[colour_func_dofmap[iC]]
            cD = colour_func_arr[colour_func_dofmap[iD]]
            
            # Volumes of the two cells
            #vC = cell_info[iC].volume
            #vD = cell_info[iD].volume
                    
            # Compressive velocity
            
            #Uc = compr_fac*norm(ump)*unity_gradient
            
            # Volume to flux
            #volDelta = numpy.dot(Uc, normal)*area
            #if volDelta < 0:
            #    # This should be just noise
            #    assert volDelta < EPS*100
            #    volDelta = 0.0
            
            unity_gradient = gradient/(norm(gradient) + EPS)
            w = abs(numpy.dot(unity_gradient, normal))*compr_fac
            
            # Take no more than what exists and do not overfill
            cDelta = cC*w #max(volDelta, cC)
            cDelta = min(cDelta, 1 - cD)
            
            if cDelta < 0:
                cDelta = max(cDelta, cC-1)
                
            # Local sharpening of the colour function in D and C cells
            colour_func_arr[colour_func_dofmap[iC]] -= cDelta
            colour_func_arr[colour_func_dofmap[iD]] += cDelta
        
        self.colour_function.vector().set_local(colour_func_arr)
         
