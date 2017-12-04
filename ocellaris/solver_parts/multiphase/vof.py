import numpy
from dolfin import Constant, FunctionSpace


class VOFMixin(object):
    """
    This is a mixin class to avoid having duplicates of the methods calculating
    rho, nu and mu. Any subclass using this mixin must define the method
    "get_colour_function(k)" and can also redefine the boolean property that
    controls the way mu is calculated, "calculate_mu_directly_from_colour_function".
    """
    calculate_mu_directly_from_colour_function = True
    default_polynomial_degree_colour = 0
    
    @classmethod
    def create_function_space(cls, simulation):
        mesh = simulation.data['mesh']
        cd = simulation.data['constrained_domain']
        Vc_name = simulation.input.get_value('multiphase_solver/function_space_colour',
                                             'Discontinuous Lagrange', 'string')
        Pc = simulation.input.get_value('multiphase_solver/polynomial_degree_colour',
                                        cls.default_polynomial_degree_colour, 'int')
        Vc = FunctionSpace(mesh, Vc_name, Pc, constrained_domain=cd)
        simulation.data['Vc'] = Vc
        simulation.ndofs += Vc.dim()
    
    def get_colour_function(self, k):
        """
        Return the colour function on timestep t^{n+k}
        """
        raise NotImplementedError('The get_colour_function method must be implemented by subclass!')
    
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
        return Constant(self.rho0) * c + Constant(self.rho1) * (1 - c)
    
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
        return Constant(self.nu0) * c + Constant(self.nu1) * (1 - c)
    
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
            mu0 = self.nu0 * self.rho0
            mu1 = self.nu1 * self.rho1
            return Constant(mu0) * c + Constant(mu1) * (1 - c)
        
        else:
            nu = self.get_laminar_kinematic_viscosity(k, c)
            rho = self.get_density(k, c)
            return nu * rho
    
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
            mu0 = self.nu0 * self.rho0
            mu1 = self.nu1 * self.rho1
            return min(mu0, mu1), max(mu0, mu1)
        else:
            c = numpy.linspace(0, 1, 1000)
            nu = self.nu0 * c + self.nu1 * (1 - c)
            rho = self.rho0 * c + self.rho1 * (1 - c)
            mu = nu * rho
            return mu.min(), mu.max()
    