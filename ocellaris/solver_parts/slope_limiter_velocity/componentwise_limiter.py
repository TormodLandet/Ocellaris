# encoding: utf8
from __future__ import division
from ocellaris.solver_parts.slope_limiter import SlopeLimiter
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase


@register_velocity_slope_limiter('Componentwise')
class ComponentwiseSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    description = 'Scalar limiting of each component'
    
    def __init__(self, simulation, vel, vel_name, use_cpp=True):
        """
        Use a standard slope limiter on each component of the velocity field 
        """
        self.additional_plot_funcs = []
        self.limiters = [] 
        dim, = vel.ufl_shape 
        
        # Create slope limiters
        for d in range(dim):
            name = '%s%d' % (vel_name, d)
            lim = SlopeLimiter(simulation, vel_name, vel[d], name)
            self.additional_plot_funcs.extend(lim.additional_plot_funcs)
            self.limiters.append(lim)
    
    def run(self):
        """
        Perform slope limiting of the velocity field
        """
        for lim in self.limiters:
            lim.run()
