import dolfin
from . import register_boundary_condition, BoundaryCondition
from ocellaris.utils import CodedExpression, OcellarisCppExpression


class OcellarisNeumannBC(object):
    def __init__(self, simulation, value, subdomain_id):
        """
        A simple storage class for Neumann conditions. This is
        used when defining the linear part of the weak forms
        """
        self.simulation = simulation
        self._value = value
        self.subdomain_id = subdomain_id
        
    def func(self):
        """
        The boundary value derivative function 
        """
        return self._value
    
    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)
    
    def __repr__(self):
        return '<OcellarisNeumannBC on subdomain %d>' % self.subdomain_id


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
        df_value = dolfin.Constant(value)
        
        # Store the boundary condition for use in the solver
        bc = OcellarisNeumannBC(self.simulation, df_value, subdomain_id)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)
        
        self.simulation.log.info('    ConstantGradient %r for %s' % (value, var_name))


@register_boundary_condition('CodedGradient')
class CodedNeumannBoundary(BoundaryCondition):
    description = 'A coded Neumann condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Neumann condition with coded value
        """
        self.simulation = simulation
        
        # Make a dolfin Expression object that runs the code string
        code = inp_dict['code']
        
        if isinstance(code, list):
            assert len(code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                description = 'coded gradient boundary condition for %s' % name
                expr = CodedExpression(simulation, code[0], description)
                self.register_neumann_condition(name, expr, subdomains, subdomain_id)
        else:
            description = 'coded gradient boundary condition for %s' % var_name
            expr = CodedExpression(simulation, code, description)
            self.register_neumann_condition(var_name, expr, subdomains, subdomain_id)
       
    def register_neumann_condition(self, var_name, expr, subdomains, subdomain_id):
        """
        Store the boundary condition for use in the solver
        """
        bc = OcellarisNeumannBC(self.simulation, expr, subdomain_id)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    Coded gradient for %s' % var_name)


@register_boundary_condition('CppCodedGradient')
class CppCodedNeumannBoundary(BoundaryCondition):
    description = 'A C++ coded Neumann boundary condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        C++ coded Neumann condition
        """
        self.simulation = simulation
        cpp_code = inp_dict['cpp_code']
        
        if isinstance(cpp_code, list):
            assert len(cpp_code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_neumann_condition(name, cpp_code[d], subdomain_id)
        else:
            self.register_neumann_condition(var_name, cpp_code, subdomain_id)
    
    def register_neumann_condition(self, var_name, cpp_code, subdomain_id):
        """
        Add a C++ coded Neumann condition to this variable
        """
        description = 'boundary condititon for %s' % var_name
        expr = OcellarisCppExpression(self.simulation, cpp_code, description, update=True)
        bc = OcellarisNeumannBC(self.simulation, expr, subdomain_id)
        bcs = self.simulation.data['neumann_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    C++ coded gradient for %s' % var_name)