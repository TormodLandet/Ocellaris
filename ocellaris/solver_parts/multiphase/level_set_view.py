import dolfin
from ocellaris.probes.free_surface_locator import get_free_surface_locator


class LevelSetView:
    def __init__(self, simulation):
        """
        A LevelSetView is a view of another multi-phase model (typically
        VOF) as a level set function. This can be handy when distance from
        the free surface is needed in the code.

        This implementation of a level set function is not a stand-alone
        multi-phase model and cannot itself advect the free surface /
        density field.
        """
        self.simulation = simulation
        self.base = None
        self.base_type = 'unset'
        self.base_name = None
        self.iso_value = None
        self.name = None

        # Create the level set function
        mesh = simulation.data['mesh']
        V = dolfin.FunctionSpace(mesh, 'CG', 1)
        self.level_set_function = dolfin.Function(V)
        self.cache = preprocess(simulation, self.level_set_function)

    def set_density_field(self, c, name='c', value=0.5, update_hook='MultiPhaseModelUpdated'):
        """
        Set the generalised density field, either the vof function or
        the direct tensity field from a variable density simulation.
        This must be done before using the view
        """
        self.base = c
        self.base_type = 'vof'
        self.base_name = name
        self.iso_value = value

        # Store the level set function in an easy to access location and make
        # sure it is saved to any restart files that are written
        self.name = 'ls_%s_%s' % (name, float_to_ident(value))
        self.level_set_function.rename(self.name, self.name)
        self.simulation.data[self.name] = self.level_set_function

        self._locator = get_free_surface_locator(self.simulation, name, c, value)
        self._locator.add_update_hook(update_hook, self.update)

    def update(self):
        """
        The underlying data has changed, update the level set function
        """
        if self.base_type == 'vof':
            self._update_from_vof()
        else:
            raise NotImplementedError(
                'Cannot compute level set function ' 'from %r base field' % self.base_type
            )

    def _update_from_vof(self):
        # This can be expensive, will involve recomputing the crossing
        # points if the density function has changed since the last access
        crossings = self._locator.crossing_points
        update_level_set_view(self.simulation, self.level_set_function, crossings, self.cache)


def float_to_ident(v):
    s = repr(v)
    for s0, s1 in [('.', '_'), ('+', 'p'), ('-', 'm')]:
        s = s.replace(s0, s1)
    return s


def preprocess(simulation, level_set_view):
    V = level_set_view.function_space()


def update_level_set_view(simulation, level_set_view, crossings, cache):
    pass
