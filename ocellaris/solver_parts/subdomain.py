"""
Subdomains for Ocellaris

A subdomain in this context is a way to define a cell-wise constant
function that varies between 0 and 1 depending on whether a cell is inside
(1.0) or outside (0.0) the subdomain. The boundaries between inside and 
outside are allowed to be smooth (values in between 0.0 and 1.0).

One reason for creating smooth subdomains is to mark cells near the free
surface with a gradual falloff to 0 in cells outside the free surface
region. This can be used for boundary conditions (interface slip) and
limiters (only limit near the free surface).
"""
import dolfin
from ocellaris.utils import OcellarisError


def make_subdomain(simulation, inp):
    name = inp.get_value('name', 'unnamed', 'string')
    domain_type = inp.get_value('type', required_type='string')

    if domain_type in SUBDOMAIN_TYPES:
        return SUBDOMAIN_TYPES[domain_type](simulation, name, inp)
    else:
        valid = ', '.join(SUBDOMAIN_TYPES.keys())
        raise OcellarisError(
            "Invalid subdomain type",
            "%s is not a valid subdomain type, specify one of %s" % (domain_type, valid),
        )


class OcellarisSubdomain:
    def __init__(self, simulation, name, inp):
        """
        A subdomain has a .function attribute which is a DG0 function with
        values between 1.0 (inside) and 0.0 (outside) to mark the 
        subdomain in a way that allows smooth transitions.
        """
        self.simulation = simulation
        self.name = name

        mesh = simulation.data['mesh']
        self.V = dolfin.FunctionSpace(mesh, 'DG', 0)
        self.function = dolfin.Function(self.V)

        func_name = 'subdomain_%s' % name
        self.function.rename(func_name, func_name)
        simulation.data[func_name] = self.function

        if inp.get_value('plot', False, 'bool'):
            simulation.io.add_extra_output_function(self.function)

        self._callbacks = []
        self._initialize(inp)

    def add_update_callback(self, func):
        self._callbacks.append(func)

    def _initialize(self, inp):
        raise NotImplementedError('A subclass should implement this')


class NearInterface(OcellarisSubdomain):
    def _initialize(self, inp):
        mpm = self.simulation.multi_phase_model
        self.level_set_view = mpm.get_level_set_view()
        self.level_set_view.add_update_callback(self.update)

        # Read input
        self.radius = inp.get_value('radius', required_type='float')

        # Form to compute the cell average distace to the free surface
        v = dolfin.TestFunction(self.V)
        ls = self.level_set_view.level_set_function
        cv = dolfin.CellVolume(self.V.mesh())
        self.dist_form = dolfin.Form(ls * v / cv * dolfin.dx)

    def update(self):
        # Compute the cell average distance to the free surface
        f = self.function
        dolfin.assemble(self.dist_form, tensor=f.vector())

        # Find the inside cells and the cells in the smoothing zone
        arr = f.vector().get_local() / self.radius
        inside = arr <= 1.0
        smooth = (arr > 1.0) & (arr <= 2.0)
        outside = arr > 2.0

        # Create the inside and smoothing zone with a smooth transition
        # from 1 to 0 when r goes from 1 to 2. The slope in both ends is 0
        arr[inside] = 1.0
        r = arr[smooth]
        arr[smooth] = 2 * r ** 3 - 9 * r ** 2 + 12 * r - 4
        arr[outside] = 0.0

        f.vector().set_local(arr)
        f.vector().apply('insert')

        for func in self._callbacks:
            func()


SUBDOMAIN_TYPES = {'NearInterface': NearInterface}
