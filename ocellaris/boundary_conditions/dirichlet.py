import dolfin
from . import register_boundary_condition, BoundaryCondition

@register_boundary_condition('ConstantValue')
class DirichletBC(BoundaryCondition):
    description = 'A prescribed constant value Dirichlet condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with constant value
        """
        self.simulation = simulation
        
        value = inp_dict['value']
        if isinstance(value, list):
            assert len(value) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_dirichlet_condition(name, value[d], subdomains, subdomain_id)
        else:
            self.register_dirichlet_condition(var_name, value, subdomains, subdomain_id)
    
    def register_dirichlet_condition(self, var_name, value, subdomains, subdomain_id):
        """
        Add a Dirichlet condition to this variable
        """
        # Get the variables function space
        func = self.simulation.data[var_name]
        V = func.function_space()
        
        # Store the boundary condition for use in the solver
        assert isinstance(value, (float, int, long))
        value = dolfin.Constant(value)
        bc = dolfin.DirichletBC(V, value, subdomains, subdomain_id)
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
