import dolfin
from . import register_boundary_condition, BoundaryConditionCreator
from ocellaris.utils import CodedExpression, OcellarisCppExpression, OcellarisError


class OcellarisDirichletBC(dolfin.DirichletBC):
    def __init__(self, simulation, V, value, subdomain_marker, subdomain_id, updater=None):
        """
        A simple storage class for Dirichlet conditions. This is
        used when defining the linear part of the weak forms and
        for normal boundary strong conditions
        """
        super(OcellarisDirichletBC, self).__init__(V, value, subdomain_marker, subdomain_id, method='geometric')
        self.simulation = simulation
        self._value = value
        self.subdomain_marker = subdomain_marker
        self.subdomain_id = subdomain_id
        self._updater = updater
    
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
    
    def copy_and_change_function_space(self, V):
        """
        Return a copy with a new function space. Used when converting from
        BCs for a segregated solver (default) to BCs for a coupled solver
        """
        return OcellarisDirichletBC(self.simulation, V, self._value,
                                    self.subdomain_marker, self.subdomain_id)
    
    def update(self):
        """
        Update the time and other parameters used in the BC.
        This is used every timestep and for all RK substeps
        """
        if self._updater:
            self._updater(self.simulation.timestep,
                          self.simulation.time,
                          self.simulation.dt)
    
    def __repr__(self):
        return '<OcellarisDirichletBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('ConstantValue')
class ConstantDirichletBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant value Dirichlet condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with constant value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]
        
        value = inp_dict.get_value('value', required_type='any')
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
        if not isinstance(value, (float, int)):
            raise OcellarisError('Error in ConstantValue BC for %s' % var_name,
                                 'The value %r is not a number' % value)
        df_value = dolfin.Constant(value)
        
        # Store the boundary condition for use in the solver
        bc = OcellarisDirichletBC(self.simulation, self.func_space, df_value, subdomains, subdomain_id)
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        
        self.simulation.log.info('    Constant value %r for %s' % (value, var_name))


@register_boundary_condition('CodedValue')
class CodedDirichletBoundary(BoundaryConditionCreator):
    description = 'A coded Dirichlet condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with coded value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]
        
        # Make a dolfin Expression object that runs the code string
        code = inp_dict.get_value('code', required_type='any')
        
        if isinstance(code, list):
            assert len(code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                description = 'coded value boundary condition for %s' % name
                sub_code = inp_dict.get_value('code/%d' % d, required_type='string')
                expr = CodedExpression(simulation, sub_code, description)
                self.register_dirichlet_condition(name, expr, subdomains, subdomain_id)
        else:
            description = 'coded value boundary condition for %s' % var_name
            expr = CodedExpression(simulation, code, description)
            self.register_dirichlet_condition(var_name, expr, subdomains, subdomain_id)
    
    def register_dirichlet_condition(self, var_name, expr, subdomains, subdomain_id):
        """
        Store the boundary condition for use in the solver
        """
        bc = OcellarisDirichletBC(self.simulation, self.func_space, expr, subdomains, subdomain_id)
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    Coded value for %s' % var_name)


@register_boundary_condition('CppCodedValue')
class CppCodedDirichletBoundary(BoundaryConditionCreator):
    description = 'A C++ coded Dirichlet condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with C++ coded value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]
        
        # Make a dolfin Expression object that runs the code string
        code = inp_dict.get_value('cpp_code', required_type='any')
        
        if isinstance(code, list):
            assert len(code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                sub_code = inp_dict.get_value('cpp_code/%d' % d, required_type='string')
                self.register_dirichlet_condition(name, sub_code, subdomains, subdomain_id)
        else:
            self.register_dirichlet_condition(var_name, code, subdomains, subdomain_id)
    
    def register_dirichlet_condition(self, var_name, cpp_code, subdomains, subdomain_id):
        """
        Store the boundary condition for use in the solver
        """
        description = 'boundary condititon for %s' % var_name
        P = self.func_space.ufl_element().degree()
        expr, updater = OcellarisCppExpression(self.simulation, cpp_code, description,
                                               P, return_updater=True)
        
        bc = OcellarisDirichletBC(self.simulation, self.func_space, expr,
                                  subdomains, subdomain_id, updater=updater)
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    C++ coded value for %s' % var_name)
