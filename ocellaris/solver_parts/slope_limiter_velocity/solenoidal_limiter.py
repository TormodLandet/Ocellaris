# encoding: utf8
from __future__ import division
from solenoidal import SolenoidalLimiter
from ocellaris.utils import verify_key
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase


@register_velocity_slope_limiter('Solenoidal')
class SolenoidalSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'Limit in a divergence free polynomial space'
    
    def __init__(self, simulation, vel, vel_name, use_cpp=True):
        """
        Use a solenoidal polynomial slope limiter on the velocity field 
        """
        # Verify input
        V = vel[0].function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        loc = 'SolenoidalSlopeLimiterVelocity'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (2,), loc)
        verify_key('function shape', vel.ufl_shape, [(2,)], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        
        # Store input
        self.simulation = simulation
        self.vel = vel
        self.vel_name = vel_name
        self.degree = degree
        self.mesh = mesh
        self.use_cpp = use_cpp
        
        # Create slope limiter
        self.sollim = SolenoidalLimiter(vel)
    
    def run(self):
        """
        Perform slope limiting of the velocity field
        """
        self.sollim.run()
