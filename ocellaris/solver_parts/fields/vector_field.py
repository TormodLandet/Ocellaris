import dolfin
from ocellaris.utils import ocellaris_error, OcellarisCppExpression
from . import register_known_field, KnownField, DEFAULT_POLYDEG


@register_known_field('VectorField')
class VectorField(KnownField):
    description = 'Vector field'
    
    def __init__(self, simulation, field_inp):
        """
        A vector field
        """
        self.simulation = simulation
        self.read_input(field_inp)
        
        # Show the input data
        simulation.log.info('Creating a vector field %r' % self.name)
        simulation.log.info('    Variable: %r' % self.var_name)
        simulation.log.info('    Stationary: %r' % self.stationary)
        simulation.log.info('    Polynomial degree: %r' % self.polydeg)
        
        self.exprs = None
        self.funcs = None
        self.updaters = None
        self.V = dolfin.FunctionSpace(simulation.data['mesh'], 'CG', self.polydeg)
        simulation.hooks.add_pre_timestep_hook(self.update, 'Update scalar field %r' % self.name)
    
    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        self.var_name = field_inp.get_value('variable_name', 'phi', required_type='string')
        self.stationary = field_inp.get_value('stationary', False, required_type='bool')
        self.polydeg = field_inp.get_value('polynomial_degree', DEFAULT_POLYDEG, required_type='int')
        self.cpp_code = field_inp.get_value('cpp_code', required_type='list(string)')
        if len(self.cpp) != self.simulation.ndim:
            ocellaris_error('Wrong number of dimensions',
                            '')
    
    def update(self, timestep_number, t, dt):
        """
        Called by simulation.hooks on the start of each time step
        """
        if self.stationary:
            return
        for u in self.updaters:
            u(timestep_number, t, dt)
    
    def _get_expression(self):
        if self.exprs is None:
            self.exprs = []
            self.updaters = []
            for d, cpp in enumerate(self.cpp_code):
                expr, updater = OcellarisCppExpression(self.simulation, cpp,
                                                       'Vector field %s component %d' % (self.name, d),
                                                       self.polydeg, update=False,
                                                       return_updater=True)
                self.exprs.append(expr)
                self.updaters.append(updater)
            self.exprs = expr
            self.updaters = updater
        return self.exprs
    
    def get_variable(self, name):
        if not name == self.var_name:
            ocellaris_error('Vector field does not define %r' % name,
                            'This vector field defines %r' % self.var_name)
        if self.funcs is None:
            funcs = []
            for e in self._get_expression():
                func = dolfin.Function(self.V)
                func.interpolate(e)
                funcs.append(func)
            self.funcs = dolfin.as_vector(funcs)
        return self.funcs
