import dolfin
from ocellaris.utils import report_error, RunnablePythonString

class BoundaryRegion(object):
    def __init__(self, simulation, marker, index):
        """
        Create boundary conditions for the given part
        
        This class reads the input for the boundary part with the given index
        and creates boundary condition objects for each of the listed functions. 
        
        Arguments:
            simulation: The simulation object
            marker: a facet function used to mark the boundary
            index: the number of this part in the list of boundary conditions
                dictionaries in the simulation input. The mark in the marker
                function will be this number plus one	
        """
        self.simulation = simulation
        self.marker = marker
        self.index = index
        self.conditions = {}
        
        inp = simulation.input.get_value('boundary_conditions', required_type='list(dict)')[index]
        self.name = inp['name']
        self.selector_name = inp['selector']
        self.input = inp
        
        # Let the default boundary marking 0 be for unmarked regions
        self.mark_id = index + 1
        
        # Mark the region of the boundary covered by this boundary condition
        if self.selector_name == 'region':
            self.selector = RegionSelector(simulation)
            code_string = inp['inside_code']
            self.selector.set_inside_code(code_string, self.name)
            try:
                self.selector.mark(marker, self.mark_id)
            except Exception as e:
                report_error('Error in boundary condition',
                             'Marking boundary "%s" with region_code="%s" failed. '
                             % (self.name, code_string) +
                             '\n\nThe error was "%s"' % e +
                             '\n\nDid you remember that x is an array?')
        else:
            report_error('Error: unknown boundary selector',
                         'Boundary condition for boundary "%s" has '
                         'selector="%s". This selector is not implemented.'
                         '\n\nImplemented selectors:\n\n'
                         ' - region'
                         % (self.name, self.selector_name))
    
    def create_periodic_boundary_conditions(self):
        """
        Create periodic boundary conditions for this region
        
        This is separated from the normal boundary conditions
        because Dirichlet boundary conditions depend on having
        function spaces created, while the function spaces themselves
        are dependent on the definition of periodic boundaries. 
        """
        sim = self.simulation
        if 'map_code' in self.input:
            sim.log.info('Applying periodic boundary conditions on %s' %  self.name)
            
            if self.simulation.data['constrained_domain'] is not None:
                report_error('Error in specification of periodic boundary conditions',
                             'There can only be one periodic boundary region in the domain. '
                             'Found more than one periodic region when processing boundary conditions')
            
            self.selector.set_map_code(self.input['map_code'], self.name)
            self.simulation.data['constrained_domain'] = self.selector
        
    def create_boundary_conditions(self):
        """
        Create the Dirichlet and Neumann boundary conditions for this region
        """
        sim = self.simulation
        
        # Get boundary conditions on this boundary
        for key, value in self.input.items():
            if not isinstance(value, dict) or 'type' not in value:
                continue
            
            bc_type = value['type']
            sim.log.info('Applying %s boundary condition for %s on %s' % 
                         (bc_type, key, self.name))
            bc_class = get_boundary_condition(bc_type)
            bc = bc_class(sim, key, value, self.marker, self.mark_id)
            self.conditions[key] = bc


class RegionSelector(dolfin.SubDomain):
    def __init__(self, simulation):
        """
        A sub domain used to mark the boundaries of the domain in order
        to apply boundary conditions. The code that runs in the inside()
        and map() methods is provided by the user on the input file
        """
        super(RegionSelector, self).__init__()
        self.simulation = simulation
        
    def set_inside_code(self, code_string, region_name):
        """
        Set the code to be used in the .inside() method to define
        the region
        
        It the code contains a newline it must be the core
        of the function and define the inside variable
        """
        self.inside_func = RunnablePythonString(self.simulation, code_string,
                                                'inside code for %s' % region_name,
                                                var_name='inside')
    
    def set_map_code(self, code_string, region_name):
        """
        Set the code to be used in the .map() method for periodic
        boundary conditions.
        
        It the code contains a newline it must be the core
        of the function and also include a return statement,
        otherwise it is assumed to be an expression
        """
        self.map_func = RunnablePythonString(self.simulation, code_string,
                                             'map code for %s' % region_name)
    
    def inside(self, x, on_boundary):
        """
        Determine whether the point x is in this region
        """
        return self.inside_func.run(x=x, on_boundary=on_boundary)
    
    def map(self, x, y):
        """
        Map from x[i] to y[i] on periodic boundary 
        """
        return self.map_func.run(x=x, y=y)

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
    
    def func(self):
        """
        Returns the value at the boundary for Dirichlet boundary conditions
        and the normal derivative at the boundaru for Neumann bcs. 
        """
        raise NotImplementedError()
    
    def ds(self):
        """
        Returns the ds measure of the part of the boundary which this boundary
        condition applies to
        """
        raise NotImplementedError()

from . import dirichlet
from . import neumann
