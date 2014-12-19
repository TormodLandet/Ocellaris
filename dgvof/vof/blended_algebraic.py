from dolfin import FunctionSpace, Function, DirichletBC, \
                   TrialFunction, TestFunction, Constant, \
                   FacetNormal, dx, ds, dS, \
                   inner, grad, lhs, rhs, jump, solve
from dgvof.convection import get_convection_scheme

class BlendedAlgebraicVofScheme():
    def __init__(self, mesh, convection_scheme, velocity_field):
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
        # Define function space and solution function
        self.function_space = FunctionSpace(mesh, "DG", 0)
        self.colour_function = Function(self.function_space)
        self.prev_colour_function = cp = Function(self.function_space)
        
        # Create test and trial functions
        c = TrialFunction(self.function_space)
        v = TestFunction(self.function_space)
        
        # The convection blending function that counteracts numerical diffusion
        scheme = get_convection_scheme(convection_scheme)
        self.convection_scheme = scheme(mesh, self.function_space)
        beta = self.convection_scheme.blending_function
        
        # The time step (real value to be supplied later)
        self.dt = Constant(1.0)
        
        # The normal on each face
        normal = FacetNormal(mesh)
        
        # Upstream and downstream normal velocities
        vel = velocity_field
        vel_nU = (inner(vel, normal) + abs(inner(vel, normal)))/2
        vel_nD = (inner(vel, normal) - abs(inner(vel, normal)))/2
        
        # The blended flux
        flux = jump((1-beta)*vel_nU*c + beta*vel_nD*c)
        
        # Equation to solve
        eq = (c-cp)/self.dt*v*dx \
             - c*inner(vel, grad(v))*dx \
             + flux*jump(v)*dS \
             + c*inner(vel, normal)*v*ds
        self.lhs = lhs(eq)
        self.rhs = rhs(eq)
        
        # Boundary condition (C=0 on boundary)
        # TODO: this is not very general ...
        self.bc = DirichletBC(self.function_space, 0, lambda x, on_boundary: on_boundary)
        
        # Store the reference to the convecting velocity field
        self.velocity_field = velocity_field
    
    def update(self, t, dt):
        """
        Update the VOF field by advecting it for a time dt
        using the given divergence free velocity field
        """
        # Update the convection blending factors
        self.convection_scheme.update(t, dt, self.colour_function, self.velocity_field)
        
        # Solve the advection equation
        self.dt.assign(dt)
        solve(self.lhs == self.rhs, self.colour_function, self.bc)
        self.prev_colour_function.assign(self.colour_function)
