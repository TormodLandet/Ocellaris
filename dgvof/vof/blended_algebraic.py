import numpy
import dolfin
from dolfin import FunctionSpace, Function, DirichletBC, \
                   TrialFunction, TestFunction, Constant, \
                   FacetNormal, dx, ds, dS, solve, \
                   inner, grad, lhs, rhs, jump
from dgvof.convection import get_convection_scheme

class BlendedAlgebraicVofScheme():
    def __init__(self, simulation):
        """
        A blended alegraic VOF scheme works by using a specific 
        convection scheme in the advection of the colour function
        that ensures a sharp interface.
        
        * The mesh should be the same as is used for the velocities
        * The convection scheme should be the name of a convection
          scheme that is tailored for advection of the colour 
          function, i.e "HRIC", "MHRIC", "RHRIC" etc, 
        * The velocity field should be divergence free
        """
        self.simulation = simulation
        self.mesh = simulation.data['mesh']
        
        # Define function space and solution function
        self.function_space = FunctionSpace(self.mesh, "DG", 0)
        self.colour_function = Function(self.function_space)
        self.prev_colour_function = cp = Function(self.function_space)
        
        simulation.data['c'] = self.colour_function
        simulation.data['c_p'] = self.prev_colour_function
        simulation.input.setdefault('convection', {}).setdefault('c', {})['compute_facet_gradient'] = True
        
        # Create test and trial functions
        c = TrialFunction(self.function_space)
        v = TestFunction(self.function_space)
        
        # The convection blending function that counteracts numerical diffusion
        scheme = get_convection_scheme(simulation.input['convection']['c'].get('convection_scheme', 'HRIC'))
        self.convection_scheme = scheme(simulation, 'c')
        beta = self.convection_scheme.blending_function
        
        # The time step (real value to be supplied later)
        self.dt = Constant(1.0)
        
        # The normal on each face
        normal = FacetNormal(self.mesh)
        
        # Upstream and downstream normal velocities
        vel = simulation.data['u']
        flux_nU = c*(inner(vel, normal) + abs(inner(vel, normal)))/2
        flux_nD = c*(inner(vel, normal) - abs(inner(vel, normal)))/2

        # Define the blended flux
        # The blending factor beta is not DG, so beta('+') == beta('-')
        b = beta('+')
        flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
        
        # Reconstruct the gradient from the colour function DG0 field
        self.convection_scheme.gradient_reconstructor.initialize()
        gradient = self.convection_scheme.gradient_reconstructor.gradient
        
        # Compressive flux
        #compr_fac = Constant(simulation.input['VOF'].get('compression_factor', 1.0))
        #norm = lambda f: (f.sub(0)**2 + f.sub(1)**2)**0.5
        #compr_vel = compr_fac * norm(vel) * gradient / (norm(gradient) + Constant(1e-6))
        #cvelU = (inner(compr_vel, normal) + abs(inner(compr_vel, normal)))/2
        #q = c * (1-cp)
        #cflux = q('+')*cvelU('+') - q('-')*cvelU('-')
        compr_vel = Constant((0.0, 0.0))
        cflux = Constant(0.0)
        
        # Equation to solve
        eq = (c-cp)/self.dt*v*dx \
             - c*inner(vel, grad(v))*dx \
             - c*(1-cp)*inner(compr_vel, grad(v))*dx \
             + (flux + cflux)*jump(v)*dS \
             + c*inner(vel, normal)*v*ds
        self.lhs = lhs(eq)
        self.rhs = rhs(eq)
        
        # Boundary condition (C=0 on boundary)
        # TODO: this is not very general ...
        self.bc = DirichletBC(self.function_space, 0, lambda x, on_boundary: on_boundary)
        
        # Store the reference to the convecting velocity field
        self.velocity_field = simulation.data['u']
        
        simulation.plotting.add_plot('c', self.colour_function)
        simulation.plotting.add_plot('c_grad', gradient)
        #import dolfin
        #simulation.add_plot('cvel', dolfin.project(cvel, V=fgrad.function_space()))
        #print 'cvel'
        #self.plots.append(('cvel', Plot2DCR1Vec()))
        #print 'cf'
        #self.plots.append(('cf', Plot2DCR1Vec(dolfin.project(cf, V=fgrad.function_space()))))
        
    def update(self, t, dt):
        """
        Update the VOF field by advecting it for a time dt
        using the given divergence free velocity field
        """
        # Reconstruct the gradient
        self.convection_scheme.gradient_reconstructor.reconstruct()
        
        # Update the convection blending factors
        self.convection_scheme.update(t, dt, self.velocity_field)
        
        # Solve the advection equation
        self.dt.assign(dt)
        solve(self.lhs == self.rhs, self.colour_function, self.bc)
        self.prev_colour_function.assign(self.colour_function)
        
        # Compress interface
        self.compress(t, dt)
    
    def compress(self, t, dt):
        """
        Explicit compression
        """
        #self.simulation.plotting.plot('c', '_uncompr')
        
        compr_fac = self.simulation.input['VOF'].get('compression_factor', 1.0)
        if compr_fac == 0:
            return
        
        facet_info = self.simulation.data['facet_info']
        #cell_info = self.simulation.data['cell_info']
        conFC = self.simulation.data['connectivity_FC']
        ndim = self.simulation.ndim
        
        colour_func_vec = self.colour_function.vector()
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
            c0 = colour_func_vec[colour_func_dofmap[i0]]
            c1 = colour_func_vec[colour_func_dofmap[i1]]
            
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
            cC = colour_func_vec[colour_func_dofmap[iC]]
            cD = colour_func_vec[colour_func_dofmap[iD]]
            
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
            colour_func_vec[colour_func_dofmap[iC]] -= cDelta
            colour_func_vec[colour_func_dofmap[iD]] += cDelta            
