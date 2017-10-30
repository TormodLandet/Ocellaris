import numpy
import dolfin
from ocellaris.utils import facet_dofmap, timeit, OcellarisCppExpression, OcellarisError
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
            self.simulation.log.info('    Found existing slip factor %r for %r. Reusing that'
                                     % (self.slip_factor_name, self.bc_var_name))
            return self.slip_factor_name
        
        if self.slip_factor_degree != 0:
            raise NotImplementedError('Slip factor must currently be piecewice constant')
        
        # Create the slip factor field
        self.scalar_field = sim.data[self.scalar_field_name]
        mesh = self.scalar_field.function_space().mesh()
        V = dolfin.FunctionSpace(mesh, 'DGT', 0)
        self.slip_factor = dolfin.Function(V)
        sim.data[self.slip_factor_name] = self.slip_factor
        
        # Get external facing facet midpoints and neighbours
        self.simulation = sim
        self.preprocess_facets()
        
        # Update the field before each time step
        sim.hooks.add_custom_hook(self.custom_hook_point, self.update,
                                  'Update slip length "%s"' % self.slip_factor_name)
        return self.slip_factor_name
    
    def preprocess_facets(self):
        sim = self.simulation
        conn_FV = sim.data['connectivity_FV']
        conn_VF = sim.data['connectivity_VF']
        
        self.facets = [(fidx, f) for fidx, f in enumerate(sim.data['facet_info']) if f.on_boundary]
        self.facet_dofs = facet_dofmap(self.slip_factor.function_space())
        
        # For each boundary facet find the neighbour external facets
        external_facets = set(fidx for fidx, _ in self.facets) 
        self.facet_neighbours = {}
        for fidx, _facet in self.facets:
            self.facet_neighbours[fidx] = nbs = []
            vs = conn_FV(fidx)
            for v in vs:
                for nb in conn_VF(v):
                    if nb != fidx and nb in external_facets:
                        nbs.append(nb)
    
    @timeit.named('SlipFactorUpdater')
    def update(self):
        phi = self.scalar_field
        fac = self.slip_factor
        ls = self.scalar_field_level_set
        D = self.slip_factor_distance
        
        # Get the value of the colur function in the midpoint of each external facet
        values = {}
        val = numpy.zeros(1, float)
        for fidx, facet in self.facets:
            mp = facet.midpoint
            phi.eval(val, mp)
            values[fidx] = val[0]
        
        # Find where the level set touches the boundary
        intersections = []
        for fidx, facet in self.facets:
            nbmin, nbmax = 1e100, -1e100
            for nb in self.facet_neighbours[fidx]:
                v = values[nb]
                nbmin = min(nbmin, v)
                nbmax = max(nbmax, v)
            
            v = values[fidx]
            if v <= ls and nbmax >= ls:
                intersections.append(facet.midpoint)
        
        # Update the slip factor
        arr = fac.vector().get_local()
        arr[:] = 0.0
        if intersections:
            for fidx, facet in self.facets:
                # Find the distance to the closest intersection
                min_dist = 1e100
                mp = facet.midpoint
                for pos in intersections:
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
        fac.vector().set_local(arr)
        fac.vector().apply('insert')
        #from matplotlib import pyplot
        #from ocellaris.utils.plotting_trace import plot_matplotlib_dgt
        #c = plot_matplotlib_dgt(fac)
        #pyplot.colorbar(c)
        #pyplot.savefig('debug.png')
