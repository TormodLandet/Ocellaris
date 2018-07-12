from ocellaris import Simulation

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

output:
    log_enabled: no
    stdout_enabled: no
    solution_properties: off
    save_restart_file_at_end: off
"""


def test_level_set_view():
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    sim.input.set_value('initial_conditions/cp/cpp_code', 'x[1] > 0.5 ? 0.0 : 1.0')
    sim.log.setup()
    sim.setup()

    ls = sim.multi_phase_model.get_level_set_view()

    # First, test the surface locator
    loc = ls._locator
    assert loc._crossing_points is None
    cp = loc.crossing_points
    assert loc._crossing_points is not None
    sim.hooks.run_custom_hook('MultiPhaseModelUpdated')
    assert loc._crossing_points is None
    print(cp)
    for points in cp.values():
        for pt in points:
            print(pt)
