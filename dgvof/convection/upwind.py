"""
The HRIC upwind/downwind blending sheme
"""
from . import ConvectionScheme, register_convection_scheme

@register_convection_scheme('Upwind')
class ConvectionSchemeUpwind(ConvectionScheme):
    def __init__(self, simulation, func_name):
        """
        Implementation of the upwind convection scheme
        """
        super(ConvectionSchemeUpwind, self).__init__(simulation, func_name)
        
        # Set downwind factor to 0.0
        self.blending_function.vector()[:] = 0.0
        
    def update(self, t, dt, velocity):
        """
        Update the values of the blending function beta at the facets
        """
        # Reconstruct the gradient to calculate upstream values
        self.gradient_reconstructor.reconstruct()
