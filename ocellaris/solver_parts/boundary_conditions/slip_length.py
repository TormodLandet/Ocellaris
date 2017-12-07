import numpy
import dolfin
from ocellaris.utils import (facet_dofmap, sync_arrays, get_local, set_local,
                             timeit, OcellarisCppExpression, OcellarisError)
from .robin import OcellarisRobinBC
from . import register_boundary_condition, BoundaryConditionCreator


def df_wrap(val, description, degree, sim):
    """
    Wrap numbers as dolfin.Constant and strings as
    dolfin.Expression C++ code. Lists must be ndim
    long and contain either numbers or strings 
    """
    if isinstance(val, (int, float)):
        # A real number
        return dolfin.Constant(val)
    elif isinstance(val, str):
        # A C++ code string
        return OcellarisCppExpression(sim, val, description, degree)
    elif isinstance(val, (list, tuple)):
        D = sim.ndim
        L = len(val)
        if L != D:
            raise OcellarisError('Invalid length of list',
                                 'BC list in "%r" must be length %d, is %d.'
                                 % (description, D, L))
        
        if all(isinstance(v, str) for v in val):
            # A list of C++ code strings
            return OcellarisCppExpression(sim, val, description, degree)
        else:
            # A mix of constants and (possibly) C++ strings
            val = [df_wrap(v, description + ' item %d' % i, degree, sim) for i, v in enumerate(val)]
            return dolfin.as_vector(val)


@register_boundary_condition('SlipLength')
class SlipLengthBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant slip length (Navier) boundary condition'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Wall slip length (Navier) boundary condition with constant value
        """
        self.simulation = simulation
        vn = var_name[:-1] if var_name[-1].isdigit() else var_name
        self.func_space = simulation.data['V%s' % vn]
        dim = self.func_space.num_sub_spaces()
        default_base = 0.0 if dim == 0 else [0.0] * dim
        
        length = inp_dict.get_value('slip_length', required_type='any')
        base = inp_dict.get_value('value', default_base, required_type='any')
        self.register_slip_length_condition(var_name, length, base, subdomains, subdomain_id)
    
    def register_slip_length_condition(self, var_name, length, base, subdomains, subdomain_id):
        """
        Add a Robin boundary condition to this variable
        """
        degree = self.func_space.ufl_element().degree()
        df_blend = df_wrap(length, 'slip length for %s' % var_name, degree, self.simulation)
        df_dval = df_wrap(base, 'boundary condition for %s' % var_name, degree, self.simulation)
        df_nval = 0.0
        
        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(self.simulation, self.func_space, df_blend, df_dval, 
                              df_nval, subdomains, subdomain_id)
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)
        
        self.simulation.log.info('    Constant slip length = %r (base %r) for %s'
                                 % (length, base, var_name))


@register_boundary_condition('InterfaceSlipLength')
class InterfaceSlipLengthBoundary(BoundaryConditionCreator):
    description = 'A variable slip length (Navier) boundary condition changing along the boundary'
    
    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Wall slip length (Navier) boundary condition where the slip length is multiplied
        by a slip factor ∈ [0, 1] that varies along the domain boundary depending on the
        distance to an interface (typically a free surface between two fluids).
        """
        self.simulation = simulation
        vn = var_name[:-1] if var_name[-1].isdigit() else var_name
        self.func_space = simulation.data['V%s' % vn]
        dim = self.func_space.num_sub_spaces()
        default_base = 0.0 if dim == 0 else [0.0] * dim
        
        # Create the slip length factor and the functionality that
        # updates this factor automatically before each time step,
        # just after the multiphase solver is done determining the
        # interface position
        sfu = SlipFactorUpdater(inp_dict, var_name)
        factor_name = sfu.register(simulation)
        
        fac = simulation.data[factor_name]
        length = inp_dict.get_value('slip_length', required_type='any')
        base = inp_dict.get_value('value', default_base, required_type='float')
        self.register_slip_length_condition(var_name, length, fac, factor_name, base, subdomains, subdomain_id)
    
    def register_slip_length_condition(self, var_name, length, factor, factor_name, base, subdomains, subdomain_id):
        """
        Add a Robin boundary condition to this variable
        """
        degree = self.func_space.ufl_element().degree()
        df_length = df_wrap(length, 'slip length for %s' % var_name, degree, self.simulation)
        df_blend = df_length * factor
        df_dval = df_wrap(base, 'boundary condition for %s' % var_name, degree, self.simulation)
        df_nval = 0.0
        
        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(self.simulation, self.func_space, df_blend, df_dval, 
                              df_nval, subdomains, subdomain_id)
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)
        
        self.simulation.log.info('    Variable slip length %r (base %r) for %s'
                                 % (factor_name, base, var_name))


class SlipFactorUpdater():
    def __init__(self, inp_dict, var_name):
        """
        This class makes sure the variable slip length is updated
        on the start of every timestep, after the multi phase flow
        density field has been updated
        """
        self.bc_var_name = var_name
        self.slip_factor_distance = inp_dict.get_value('slip_factor_distance', required_type='float')
        self.slip_factor_degree = inp_dict.get_value('slip_factor_degree', 0, required_type='int')
        self.slip_factor_name = inp_dict.get_value('slip_factor_name', 'slip_factor', required_type='string')
        self.scalar_field_level_set = inp_dict.get_value('scalar_field_level_set', 0.5, required_type='string')
        self.scalar_field_name = inp_dict.get_value('scalar_field', 'c', required_type='string')
        self.custom_hook_point = inp_dict.get_value('custom_hook', 'MultiPhaseModelUpdated', required_type='string')
    
    def register(self, sim):
        if self.slip_factor_name in sim.data:
            sim.log.info('    Found existing slip factor %r for %r. Reusing that'
                         % (self.slip_factor_name, self.bc_var_name))
            return self.slip_factor_name
        
        if self.slip_factor_degree != 0:
            raise NotImplementedError('Slip factor must currently be piecewice constant')
        
        # Create the slip factor field
        scalar_field = sim.data[self.scalar_field_name]
        mesh = scalar_field.function_space().mesh()
        V = dolfin.FunctionSpace(mesh, 'DGT', 0)
        self.facet_dofs = facet_dofmap(V)
        sim.data[self.slip_factor_name] = self.slip_factor = dolfin.Function(V)
        
        # We need to find where the interface intersects the boundary
        self.intersector = BoundaryLevelSetIntersector(sim, scalar_field, self.scalar_field_level_set)
        
        # Update the field before each time step
        sim.hooks.add_custom_hook(self.custom_hook_point, self.update,
                                  'Update slip length "%s"' % self.slip_factor_name)
        return self.slip_factor_name
    
    @timeit.named('SlipFactorUpdater')
    def update(self):
        fac = self.slip_factor
        D = self.slip_factor_distance
        
        intersection_points = self.intersector.get()
        
        # Initialize the factor to 0 (far away from the interface)
        arr = get_local(fac)
        arr[:] = 0.0
        if not intersection_points:
            set_local(fac, arr, apply='insert')
            return
        
        # Update the slip factor for facets close to the interface
        for fidx, facet in self.intersector.boundary_facets:
            # Find the distance to the closest intersection by linear search and
            # not using fancier nearest neighbour techniques. Assuming that there
            # are very few intersections compared to the total number of facets 
            # this is should be much faster (not tested)
            min_dist = 1e100
            mp = facet.midpoint
            for pos in intersection_points:
                d1 = mp - pos
                d = numpy.dot(d1, d1)
                min_dist = min(min_dist, d)
            min_dist = min_dist**0.5
            
            # Update the slip factor for this facet
            dof = self.facet_dofs[fidx]
            r = min_dist/D
            if r < 1:
                arr[dof] = 1
            elif r < 2:
                # Smooth transition from 1 to 0 when r goes from 1 to 2
                # The slope in both ends is 0
                arr[dof] = 2*r**3 - 9*r**2 + 12*r - 4
            else:
                arr[dof] = 0
        
        set_local(fac, arr, apply='insert')
        #from matplotlib import pyplot
        #from ocellaris.utils.plotting_trace import plot_matplotlib_dgt
        #c = plot_matplotlib_dgt(fac)
        #pyplot.colorbar(c)
        #pyplot.savefig('debug.png')


class BoundaryLevelSetIntersector:
    def __init__(self, sim, scalar_field, level_set):
        """
        Helper class to locate where the level set (a real number)
        of a scalar field intersects the domain boundary
        """
        self.scalar_field = scalar_field
        self.level_set = level_set
        
        conn_FV = sim.data['connectivity_FV']
        conn_VF = sim.data['connectivity_VF']
        conn_FC = sim.data['connectivity_FC']
        
        # Get external facing facets
        self.boundary_facets = [(fidx, f) for fidx, f in enumerate(sim.data['facet_info'])
                                if f.on_boundary]
        
        # Scalar field info
        scalar_deg = scalar_field.ufl_element().degree()
        scalar_dm = scalar_field.function_space().dofmap()
        im = self.scalar_field.function_space().dofmap().index_map()
        num_dofs_owned = im.size(im.MapSize.OWNED)
        
        # For each boundary facet find the neighbour external facets
        # and the connected cell dof closest to the facet
        external_facets = set(fidx for fidx, _ in self.boundary_facets) 
        self.facet_neighbours = {}
        self.facet_midpoints = {}
        self.facet_scalar_dofs = {}
        self.facet_is_owned = {}
        for fidx, facet in self.boundary_facets:
            # Connected external facets
            self.facet_neighbours[fidx] = nbs = []
            vs = conn_FV(fidx)
            for v in vs:
                for nb in conn_VF(v):
                    if nb != fidx and nb in external_facets:
                        nbs.append(nb)
            self.facet_midpoints[fidx] = facet.midpoint
            
            # Connected cell dofs
            cells = conn_FC(fidx)
            assert len(cells) == 1
            cidx = cells[0]
            sdofs = scalar_dm.cell_dofs(cidx)
            if scalar_deg == 0:
                sdof = sdofs[0]
            else:
                raise NotImplementedError('Scalar field dof finder for degree '
                                          '%d is not implemed yet' % scalar_deg)
            self.facet_scalar_dofs[fidx] = sdof
            self.facet_is_owned[fidx] = sdof < num_dofs_owned
    
    def get(self):
        """
        Get the current boundary intersections
        """
        phi_arr = get_local(self.scalar_field)
        ls = self.level_set
        
        # Get the value of the colour function in the midpoint of each external facet
        values = {}
        for fidx, _facet in self.boundary_facets:
            sdof = self.facet_scalar_dofs[fidx]
            values[fidx] = phi_arr[sdof]
        
        # Find where the level set touches the boundary
        intersections = []
        for fidx, facet in self.boundary_facets:
            # We only consider facets we own
            if not self.facet_is_owned[fidx]:
                continue
            
            # We only check the facets that are "below" the interface
            fval = values[fidx]
            if fval > ls:
                continue
            
            # Find the maximum scalar value among all neighbour facets
            nbmax = -1e100
            nbmax_idx = None
            for nb in self.facet_neighbours[fidx]:
                v = values[nb]
                if v > nbmax:
                    nbmax = v
                    nbmax_idx = nb
            assert nbmax_idx is not None
            
            # Determine the position of the interface if a neighbour is "above" it
            # (above in value, not necessarily in spatial coordinate)
            if nbmax >= ls:
                nb_mp = self.facet_midpoints[nbmax_idx]
                f = 0
                if v != nbmax:
                    f = (ls - fval)/(nbmax - fval)
                pt = nb_mp * f + facet.midpoint * (1 - f)
                intersections.append(pt)
                
                print('Intersection at fidx', fidx)
                print('mp', facet.midpoint)
                print('nb_mp', nb_mp)
                print('fval', fval)
                print('nbmax', nbmax)
                print('nbmax_idx', nbmax_idx)
        
        # Sync all processes
        sync_arrays(intersections, sync_all=True)
        
        return intersections
