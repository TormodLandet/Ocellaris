import dolfin
from .robin import OcellarisRobinBC
from . import register_boundary_condition, BoundaryConditionCreator


@register_boundary_condition('ConstantSlipLength')
class ConstantSlipLengthBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant slip length (Navier) boundary condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Wall slip length (Navier) boundary condition with constant value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]
        
        length = inp_dict.get_value('slip_length', required_type='float')
        base = inp_dict.get_value('value', 0.0, required_type='float')
        self.register_slip_length_condition(var_name, length, base, subdomains, subdomain_id)
    
    def register_slip_length_condition(self, var_name, length, base, subdomains, subdomain_id):
        """
        Add a Dirichlet condition to this variable
        """
        df_blend = dolfin.Constant(length)
        df_dval = dolfin.Constant(base)
        df_nval = 0.0
        
        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(self.simulation, self.func_space, df_blend, df_dval, 
                              df_nval, subdomains, subdomain_id)
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)
        
        self.simulation.log.info('    Constant slip length = %r (base %r) for %s'
                                 % (length, base, var_name))
