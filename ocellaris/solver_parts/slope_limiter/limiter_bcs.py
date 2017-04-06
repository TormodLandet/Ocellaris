import numpy
from dolfin import Timer, Constant


class SlopeLimiterBoundaryConditions(object):
    BC_TYPE_NOT_ON_BOUNDARY = 0
    BC_TYPE_DIRICHLET = 1
    BC_TYPE_NEUMANN = 2
    
    def __init__(self, simulation, field_name, dof_region_marks, V):
        """
        This class helps slope limiting of cells adjacent to the boundary by
        providing either values of the function and its derivatives at the
        boundary 
        
        @param dict dof_region_marks: As created by get_dof_region_marks() 
        """
        self.simulation = simulation
        self.field_name = field_name
        self.function_space = V
        self.dim = V.dim()
        self.active = False
        self.set_dof_region_marks(dof_region_marks)
        self._allready_warned = set()
    
    def set_dof_region_marks(self, dof_region_marks):
        """
        In case boundary regions move around the marks can be updated here
        
        @param dict dof_region_marks: map from dof to list of regions
            containing the dof (multiple regions happens in corners etc) 
        """
        same_loc_dofs = get_same_loc_dofs(self.function_space)
        self.dof_region_marks = {}  # Map from dof to ONE region
        self.region_dofs = {}       # Map from region number to list of dofs
        
        for dof, regions in dof_region_marks.items():
            region = regions[-1]  # in case of multiple regions pick the last
            self.dof_region_marks[dof] = region
            self.region_dofs.setdefault(region, []).append(dof)
            
            # Treat all dofs in the same location in the same way
            for dof2 in same_loc_dofs[dof]:
                if dof2 in dof_region_marks:
                    continue
                self.dof_region_marks[dof2] = region
                self.region_dofs[region].append(dof2)
    
    def activate(self, active=True):
        self.active = active
        
    def _warn(self, message):
        if not message in self._allready_warned:
            self.simulation.log.warning(message)
            self._allready_warned.add(message)
    
    def get_bcs(self):
        """
        Get bc type and bc value for each dof in the function space
        """
        sim = self.simulation
        boundary_dof_type = numpy.zeros(self.dim, numpy.intc)
        boundary_dof_value = numpy.zeros(self.dim, float)
        
        if not self.active:
            return boundary_dof_type, boundary_dof_value
        
        # This is potentially slow, so we time this code
        timer = Timer("Ocellaris get slope limiter boundary conditions")
        
        # Collect Dirichlet BCs for this field
        dirichlet = {}
        for bc in sim.data['dirichlet_bcs'].get(self.field_name, []):
            region_number = bc.subdomain_id - 1
            dirichlet[region_number] = bc
        
        # Collect Neumann BCs for this field
        neumann = {}
        for bc in sim.data['neumann_bcs'].get(self.field_name, []):
            region_number = bc.subdomain_id - 1
            neumann[region_number] = bc
        
        regions = sim.data['boundary']
        for region_number, dofs in self.region_dofs.items():
            boundary_region = regions[region_number]
            
            # Get the BC object
            if region_number in dirichlet:
                bc_type = self.BC_TYPE_DIRICHLET
                value = dirichlet[region_number].func()
            elif region_number in neumann:
                bc_type = self.BC_TYPE_NEUMANN
                value = neumann[region_number].func()
            else:
                self._warn('WARNING: Field %s has no BC in region %s' %
                           (self.field_name, boundary_region.name))
                continue
            
            if isinstance(value, Constant):
                # Get constant value
                val = value.values()
                assert val.size == 1
                val = val[0]
                
                for dof in dofs:
                    boundary_dof_type[dof] = bc_type
                    boundary_dof_value[dof] = val
            else:
                self._warn('WARNING: Field %s has unsupported BC %r in region %s' %
                           (self.field_name, type(value), boundary_region.name))
        
        timer.stop()
        return boundary_dof_type, boundary_dof_value


def get_same_loc_dofs(V):
    """
    Return a dictionary mapping dof number to other dofs at the same
    location in space. V should obviously be a discontinuous space,
    otherwise there will not be multiple dofs in the same location
    """
    gdim = V.mesh().geometry().dim()
    dof_coordinates = V.tabulate_dof_coordinates().reshape((-1, gdim))
    
    # Map dof coordinate to dofs, this is for DG so multiple dofs
    # will share the same location
    coord_to_dofs = {}
    max_neighbours = 0
    for dof in xrange(len(dof_coordinates)):
        coord = tuple(round(x, 5) for x in dof_coordinates[dof])
        dofs = coord_to_dofs.setdefault(coord, [])
        dofs.append(dof)
        max_neighbours = max(max_neighbours, len(dofs)-1)
    
    # Loop through dofs at same coordinate and map them to each other
    same_loc_dofs = {}
    for dofs in coord_to_dofs.values():
        for dof in dofs:
            same_loc_dofs[dof] = tuple(d for d in dofs if d != dof)
    
    return same_loc_dofs
