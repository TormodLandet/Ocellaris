import dolfin
import numpy
from ocellaris import Simulation, setup_simulation
import pytest
from helpers import skip_in_parallel


ISO_INPUT = """
ocellaris:
    type: input
    version: 1.0
mesh:
    type: Rectangle
    Nx: 4
    Ny: 4
probes:
-   name: free_surface
    enabled: yes
    type: IsoSurface
    value: 0.5
    field: c
    custom_hook: MultiPhaseModelUpdated
multiphase_solver:
    type: BlendedAlgebraicVOF
    function_space_colour: DG
    polynomial_degree_colour: 0
solver: {type: AnalyticalSolution}
boundary_conditions: [{'name': 'all', 'selector': 'code', 'inside_code': 'on_boundary'}]
time: {dt: 1.0}
physical_properties: {nu0: 1.0, nu1: 1, rho0: 1, rho1: 1}
output: {log_enabled: no}
"""


# Fails in parallel with degree = 1 on CircleCI, but not locally. This code is
# only for 2D, so parallel support is not super important. Disabling for now
# TODO: figure out why this hangs on CircleCI
@skip_in_parallel
@pytest.mark.parametrize("degree", [0, 1, 2])
def test_isoline_horizontal(degree):
    sim = Simulation()
    sim.input.read_yaml(yaml_string=ISO_INPUT)
    sim.input.set_value('multiphase_solver/polynomial_degree_colour', degree)
    setup_simulation(sim)
    probe = sim.probes['free_surface']

    # Initial value with sharp interface at x[1] == 0.5
    Vc = sim.data['Vc']
    c = sim.data['c']
    dm = Vc.dofmap()
    arr = c.vector().get_local()
    for cell in dolfin.cells(sim.data['mesh']):
        cell_value = 1 if cell.midpoint().y() < 0.5 else 0
        for dof in dm.cell_dofs(cell.index()):
            arr[dof] = cell_value
    c.vector().set_local(arr)
    c.vector().apply('insert')

    lines = probe.run(force_active=True)
    print('\nDegree:', degree, 'Vcdim:', Vc.dim())
    print(probe.name, probe.field_name, probe.value)
    print(len(lines))

    if sim.ncpu > 1:
        raise pytest.skip()

    for x, y in lines:
        print('x', x, '\ny', y)
        assert all(abs(y - 0.5) < 1e-12)

        # Results should be in sorted order
        xdx = numpy.diff(x)
        assert all(xdx > 0) or all(xdx < 0)

    assert len(lines) == 1


@pytest.mark.parametrize("degree", [1])
def test_isoline_circle(degree):
    sim = Simulation()
    sim.input.read_yaml(yaml_string=ISO_INPUT)
    sim.input.set_value('multiphase_solver/polynomial_degree_colour', degree)
    sim.input.set_value('mesh/Nx', 10)
    sim.input.set_value('mesh/Ny', 10)
    sim.input.set_value('initial_conditions/cp/cpp_code',
                        '1.1*pow(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2), 0.5)')
    setup_simulation(sim)

    sim.data['c'].assign(sim.data['cp'])
    probe = sim.probes['free_surface']
    lines = probe.run(force_active=True)

    if False:
        from matplotlib import pyplot
        c = dolfin.plot(sim.data['c'])
        pyplot.colorbar(c)
        for x, y in lines:
            pyplot.plot(x, y)
        pyplot.savefig('test_isoline_circle_%d.png' % degree)
        pyplot.close()

    print(probe.name, probe.field_name, probe.value)
    print(len(lines))
    for x, y in lines:
        # Check that the radius is constant
        r = ((x - 0.5)**2 + (y - 0.5)**2)**0.5
        print('x', x)
        print('y', y)
        print('dr', r - 0.5 / 1.1)
        assert all(abs(r - 0.5 / 1.1) < 5e-3)

        # Check that the line is clockwise or counter clockwise
        # for all segments, no going back and forth
        theta = numpy.arctan2(y - 0.5, x - 0.5) * 180 / numpy.pi
        theta[theta < 0] += 360
        tdt = numpy.diff(theta)
        tdt2 = tdt[abs(tdt) < 340]
        print('dt', tdt)
        assert all(tdt2 > 0) or all(tdt2 < 0)

    if sim.ncpu == 1:
        # The iso surface code is not written for full parallel support
        assert len(lines) == 1
        assert x[0] == x[-1] and y[0] == y[-1], "The loop should be closed"
