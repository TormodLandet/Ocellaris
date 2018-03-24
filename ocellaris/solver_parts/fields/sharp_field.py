import dolfin
from ocellaris.utils import ocellaris_error
from . import register_known_field, KnownField


@register_known_field('SharpField')
class SharpField(KnownField):
    description = 'A field with different values above and below a plane'
    
    def __init__(self, simulation, field_inp):
        """
        A scalar field
        """
        self.simulation = simulation
        self.read_input(field_inp)
        
        # Show the input data
        simulation.log.info('Creating a sharp field %r' % self.name)
        simulation.log.info('    Variable: %r' % self.var_name)
        
        self.V = dolfin.FunctionSpace(simulation.data['mesh'], 'DG', 0)
        self.func = dolfin.Function(self.V)
        
        # Initialise the sharp static field
        dm = self.V.dofmap()
        arr = self.func.vector().get_local()
        above, below = self.val_above, self.val_below
        xpos, ypos, zpos = self.xpos, self.ypos, self.zpos
        for cell in dolfin.cells(simulation.data['mesh']):
            mp = cell.midpoint()[:]
            dof, = dm.cell_dofs(cell.index())
            if (mp[0] < xpos and mp[1] < ypos and
                (simulation.ndim == 2 or mp[2] < zpos)):
                arr[dof] = below
            else:
                arr[dof] = above
        self.func.vector().set_local(arr)
        self.func.vector().apply('insert')
    
    def read_input(self, field_inp):
        self.name = field_inp.get_value('name', required_type='string')
        self.var_name = field_inp.get_value('variable_name', 'phi', required_type='string')
        self.val_above = field_inp.get_value('value_above', required_type='float')
        self.val_below = field_inp.get_value('value_below', required_type='float')
        self.xpos = field_inp.get_value('x', 1e100, required_type='float')
        self.ypos = field_inp.get_value('y', 1e100, required_type='float')
        self.zpos = field_inp.get_value('z', 1e100, required_type='float')
    
    def get_variable(self, name):
        if not name == self.var_name:
            ocellaris_error('Sharp field does not define %r' % name,
                            'This sharp field defines %r' % self.var_name)
        return self.func
