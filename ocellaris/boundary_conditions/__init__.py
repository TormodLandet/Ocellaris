import dolfin
from ocellaris.utils import report_error

class BoundaryRegion(object):
    def __init__(self, simulation, marker, index):
        """
        Create boundary conditions for the given part
        
        The index is the number of this part in the list of boundary
        conditions dictionaries in the simulation input
        """
        self.simulation = simulation
        self.index = index
        self.conditions = {}
        
        inp = simulation.input['boundary_conditions'][index]
        self.name = inp['name']
        self.selector_name = inp['selector']
        
        # Let the default boundary marking 0 be for unmarked regions
        mark_id = index + 1
        
        # Mark the region of the boundary covered by this boundary condition
        if self.selector_name == 'region':
            self.selector = RegionSelector()
            self.selector.set_code(inp['region_code'], index)
            try:
                self.selector.mark(marker, mark_id)
            except Exception as e:
                report_error('Error in boundary condition',
                             'Marking boundary "%s" with region_code="%s" failed. '
                             % (self.name, inp['region_code']) +
                             '\n\nThe error was "%s"' % e +
                             '\n\nDid you remember that x is an array?')
        else:
            report_error('Error: unknown boundary selector',
                         'Boundary condition for boundary "%s" has '
                         'selector="%s". This selector is not implemented.'
                         '\n\nImplemented selectors:\n\n'
                         ' - region'
                         % (self.name, self.selector_name))
    
        # Get boundary conditions on this boundary
        for key, value in inp.items():
            if not isinstance(value, dict) or 'type' not in value:
                continue
            bc_type = value['type']
            simulation.log.info('Applying %s boundary condition for %s on %s' % 
                                (bc_type, key, self.name))
            bc_class = get_boundary_condition(bc_type)
            bc = bc_class(simulation, key, value, marker, mark_id)
            self.conditions[key] = bc


class RegionSelector(dolfin.SubDomain):
    def set_code(self, code, index):
        filename = '<input-file::boundary_conditions::%d::region_code>' % index
        self.code = compile(code, filename, 'eval')
    
    def inside(self, x, on_boundary):
        return on_boundary and eval(self.code)


_BOUNDARY_CONDITIONS = {}

def add_boundary_condition(name, boundary_condition_class):
    """
    Register a boundary condition
    """
    _BOUNDARY_CONDITIONS[name] = boundary_condition_class

def register_boundary_condition(name):
    """
    A class decorator to register boundary conditions
    """
    def register(boundary_condition_class):
        add_boundary_condition(name, boundary_condition_class)
        return boundary_condition_class
    return register

def get_boundary_condition(name):
    """
    Return a boundary condition by name
    """
    try:
        return _BOUNDARY_CONDITIONS[name]
    except KeyError:
        report_error('Boundary condition "%s" not found' % name,
                     'Available boundary conditions:\n' +
                     '\n'.join('  %-20s - %s' % (n, s.description) 
                               for n, s in sorted(_BOUNDARY_CONDITIONS.items())))
        raise
    
class BoundaryCondition(object):
    description = 'No description available'

from . import dirichlet
from . import neumann
