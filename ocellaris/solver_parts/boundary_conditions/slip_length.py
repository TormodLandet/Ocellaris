import numpy
import dolfin
from ocellaris.utils import (
    facet_dofmap,
    get_local,
    set_local,
    timeit,
    OcellarisCppExpression,
    OcellarisError,
)
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
            raise OcellarisError(
                'Invalid length of list',
                'BC list in "%r" must be length %d, is %d.' % (description, D, L),
            )

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
        bc = OcellarisRobinBC(
            self.simulation, self.func_space, df_blend, df_dval, df_nval, subdomains, subdomain_id
        )
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info(
            '    Constant slip length = %r (base %r) for %s' % (length, base, var_name)
        )


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
        sfu = SlipFactorUpdater(simulation, inp_dict, var_name)
        factor_name = sfu.register(simulation)

        fac = simulation.data[factor_name]
        length = inp_dict.get_value('slip_length', required_type='any')
        base = inp_dict.get_value('value', default_base, required_type='float')
        self.register_slip_length_condition(
            var_name, length, fac, factor_name, base, subdomains, subdomain_id
        )

    def register_slip_length_condition(
        self, var_name, length, factor, factor_name, base, subdomains, subdomain_id
    ):
        """
        Add a Robin boundary condition to this variable
        """
        degree = self.func_space.ufl_element().degree()
        df_length = df_wrap(length, 'slip length for %s' % var_name, degree, self.simulation)
        df_blend = df_length * factor
        df_dval = df_wrap(base, 'boundary condition for %s' % var_name, degree, self.simulation)
        df_nval = 0.0

        # Store the boundary condition for use in the solver
        bc = OcellarisRobinBC(
            self.simulation, self.func_space, df_blend, df_dval, df_nval, subdomains, subdomain_id
        )
        bcs = self.simulation.data['robin_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info(
            '    Variable slip length %r (base %r) for %s' % (factor_name, base, var_name)
        )


class SlipFactorUpdater:
    def __init__(self, simulation, inp_dict, var_name):
        """
        This class makes sure the variable slip length is updated
        on the start of every timestep, after the multi phase flow
        density field has been updated
        """
        self.simulation = simulation
        self.bc_var_name = var_name
        self.slip_factor_distance = inp_dict.get_value(
            'slip_factor_distance', required_type='float'
        )
        self.slip_factor_degree = inp_dict.get_value('slip_factor_degree', 0, required_type='int')
        self.slip_factor_name = inp_dict.get_value(
            'slip_factor_name', 'slip_factor', required_type='string'
        )

    def register(self, sim):
        if self.slip_factor_name in sim.data:
            sim.log.info(
                '    Found existing slip factor %r for %r. Reusing that'
                % (self.slip_factor_name, self.bc_var_name)
            )
            return self.slip_factor_name

        # Create the slip factor field
        mesh = self.simulation.data['mesh']
        V = dolfin.FunctionSpace(mesh, 'DGT', 0)
        self.facet_dofs = facet_dofmap(V)
        sim.data[self.slip_factor_name] = self.slip_factor = dolfin.Function(V)

        # Update the field before each time step
        self.level_set_view = self.simulation.multi_phase_model.get_level_set_view()
        self.level_set_view.add_update_callback(self.update)

        # Form to compute the facet average distace to the free surface
        v = dolfin.TestFunction(V)
        ls = self.level_set_view.level_set_function
        fa = dolfin.FacetArea(mesh)
        self.dist_form = dolfin.Form(ls * v / fa * dolfin.ds)

        return self.slip_factor_name

    @timeit.named('SlipFactorUpdater')
    def update(self):
        fac = self.slip_factor
        D = self.slip_factor_distance

        # Initialize the factor to 0 (far away from the interface)
        arr = get_local(fac)
        arr[:] = 0.0

        # Get the level set distance at the facets
        distances = dolfin.assemble(self.dist_form)

        # Update the slip factor for facets close to the interface
        for fdof in self.facet_dofs:
            # Find the distance to the closest intersection
            min_dist = distances[fdof]

            # Update the slip factor for this facet
            r = min_dist / D
            if r < 1:
                arr[fdof] = 1
            elif r < 2:
                # Smooth transition from 1 to 0 when r goes from 1 to 2
                # The slope in both ends is 0
                arr[fdof] = 2 * r ** 3 - 9 * r ** 2 + 12 * r - 4
            else:
                arr[fdof] = 0

        set_local(fac, arr, apply='insert')
        # from matplotlib import pyplot
        # from ocellaris.utils.plotting_trace import plot_matplotlib_dgt
        # c = plot_matplotlib_dgt(fac)
        # pyplot.colorbar(c)
        # pyplot.savefig('debug.png')
