from . import register_boundary_condition, BoundaryCondition
import dolfin

class OcellarisNeumannBC(object):
    def __init__(self, simulation, value, subdomain_id):
        """
        A simple storage class for Neumann conditions. This is
        used when defining the linear part of the weak forms
        """
        self.simulation = simulation
        self.value = value
        self.subdomain_id = subdomain_id
        
    @property
    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)

@register_boundary_condition('ConstantGradient')
class NeumannBoundary(BoundaryCondition):
    description = 'A prescribed constant value Neumann condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Neumann condition
        """
        self.simulation = simulation
        value = inp_dict['value']
        
        if isinstance(value, list):
            assert len(value) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_neumann_condition(name, value[d], subdomain_id)
        else:
            self.register_neumann_condition(var_name, value, subdomain_id)
    
    def register_neumann_condition(self, var_name, value, subdomain_id):
        """
        Add a Neumann condition to this variable
        """
        assert isinstance(value, (float, int, long))
        value = dolfin.Constant(value)
        
        # Store the boundary condition for use in the solver
        bc = OcellarisNeumannBC(self.simulation, value, subdomain_id)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)
