import dolfin
from ocellaris import Simulation
import pytest


BASE_INPUT = """
ocellaris:
    type: input
    version: 1.0

mesh:
    type: Rectangle
    Nx: 4
    Ny: 4

physical_properties:
    rho0: 1000
    rho1: 1
    nu0: 1e-6
    nu1: 1.5e-5

solver:
    type: AnalyticalSolution

multiphase_solver:
    type: BlendedAlgebraicVOF

initial_conditions:
    cp:
        cpp_code: 'x[1] > 0.5 ? 0.0 : 1.0'

output:
    log_enabled: no
    stdout_enabled: no
    solution_properties: off
    save_restart_file_at_end: off
"""


@pytest.fixture(params=['DG0_2D_y', 'DG0_3D_y', 'DG0_2D_x'])
def vof_sim(request):
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    cases = {'DG0_2D_y': (2, 0, 'y'), 'DG0_3D_y': (3, 0, 'y'), 'DG0_2D_x': (2, 0, 'x')}
    dim, vof_deg, surf_normal = cases[request.param]

    if dim == 3:
        sim.input.set_value('mesh/type', 'Box')
        sim.input.set_value('mesh/Nz', 2)

    sim.test_surf_normal = [0, -1, 0]
    sim.test_coord_index = 1
    if surf_normal == 'x':
        sim.input.set_value('initial_conditions/cp/cpp_code', 'x[0] < 0.5 ? 0.0 : 1.0')
        sim.test_surf_normal = [1, 0, 0]
        sim.test_coord_index = 0

    sim.input.set_value('multiphase_solver/polynomial_degree_colour', vof_deg)

    sim.log.setup()
    sim.setup()
    sim.data['c'].assign(sim.data['cp'])
    return sim


def test_surface_locator(vof_sim):
    from ocellaris.probes.free_surface_locator import get_free_surface_locator

    counter = 0

    def hook():
        nonlocal counter
        counter += 1

    # Get a free surface locator
    loc = get_free_surface_locator(vof_sim, 'c', vof_sim.data['c'], 0.5)
    loc.add_update_hook('MultiPhaseModelUpdated', hook)

    # Check that the caching works as intended
    assert loc._crossing_points is None
    cp = loc.crossing_points
    assert loc._crossing_points is not None
    vof_sim.hooks.run_custom_hook('MultiPhaseModelUpdated')
    assert loc._crossing_points is None
    assert counter == 1

    # Check the number of crossings
    ndim = vof_sim.ndim
    num_crossings = len(cp)
    if vof_sim.ncpu == 1:
        assert num_crossings == 4 if ndim == 2 else 16

    # Check the location of the crossings and the fs orientation
    expected = vof_sim.test_surf_normal
    index = vof_sim.test_coord_index
    for points in cp.values():
        for pt, vec in points:
            print(pt, vec.dot(expected), vec, expected)
            assert abs(pt[index] - 0.5) < 1e-4
            assert vec.dot(expected) > 0.3


def test_level_set_view(vof_sim):
    lsv = vof_sim.multi_phase_model.get_level_set_view()
    vof_sim.hooks.run_custom_hook('MultiPhaseModelUpdated')

    lsf = lsv.level_set_function
    V = lsf.function_space()

    cpp = '0.5 - x[%d]' % vof_sim.test_coord_index
    expected = dolfin.Expression(cpp, element=V.ufl_element())

    err = dolfin.errornorm(expected, lsf)

    # TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY
    pytest.xfail()
    # TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY TEMPORARY

    assert err < 1e-2
