import dolfin
from . import register_boundary_condition, BoundaryConditionCreator


DEFAULT_BLEND = 1.0
DEFAULT_DVAL = 0.0
DEFAULT_NVAL = 0.0


class OcellarisRobinBC():
    def __init__(self, simulation, V, blend, dirichlet_value, neumann_value,
                 subdomain_marker, subdomain_id, updater=None):
        """
        A storage class for a Robin boundary conditions on the form
        
          n⋅∇φ = 1/b (φ0 - φ) + g
        
        Where b is a blending parameter and u0 and g are the Dirichlet
        and Neumann values for b → 0 and b → ∞ respectively.
        """
        self.simulation = simulation
        self._blend = blend
        self._dfunc = dirichlet_value
        self._nfunc = neumann_value
        self.subdomain_marker = subdomain_marker
        self.subdomain_id = subdomain_id
        self._updater = updater
    
    def blend(self):
        return self._blend
    
    def dfunc(self):
        return self._dfunc
    
    def nfunc(self):
        return self._nfunc
    
    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)
    
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
        return '<OcellarisRobinBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('ConstantRobin')
class ConstantRobinBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant value Robin boundary condition'
    
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
        
        blend = inp_dict.get_value('blend', DEFAULT_BLEND, required_type='float')
        dval = inp_dict.get_value('dval', DEFAULT_DVAL, required_type='float')
        nval = inp_dict.get_value('nval', DEFAULT_NVAL, required_type='float')
        if isinstance(blend, list):
            assert len(blend) == simulation.ndim
            assert len(dval) == simulation.ndim
            assert len(nval) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_robin_condition(name, blend[d], dval[d], nval[d], subdomains, subdomain_id)
        else:
            self.register_robin_condition(var_name, blend, dval, nval, subdomains, subdomain_id)
    
    def register_robin_condition(self, var_name, scale, base, grad_scale, subdomains, subdomain_id):
        """
        Add a Dirichlet condition to this variable
        """
        assert isinstance(scale, (float, int))
        assert isinstance(base, (float, int))
        assert isinstance(grad_scale, (float, int))
        df_scale = dolfin.Constant(scale)
        df_base = dolfin.Constant(base)
        df_grad_scale = dolfin.Constant(grad_scale)
        
        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(self.simulation, self.func_space, df_scale, df_base, df_grad_scale, subdomains, subdomain_id)
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)
        
        self.simulation.log.info('    Constant Robin BC %r X = %r + %r dX/dn for X = %s' % (scale, base, grad_scale, var_name))

