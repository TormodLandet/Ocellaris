# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
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
    stdout_on_all_ranks: no
    solution_properties: off
    save_restart_file_at_end: off
"""


def mk_vof_sim(case_name, modifier=None):
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    cases = {'DG0_2D_y': (2, 0, 'y'), 'DG0_3D_z': (3, 0, 'z'), 'DG0_2D_x': (2, 0, 'x')}
    dim, vof_deg, surf_normal = cases[case_name]

    if dim == 3:
        sim.input.set_value('mesh/type', 'Box')
        sim.input.set_value('mesh/Nz', 2)

    sim.test_name = case_name
    if surf_normal == 'x':
        sim.input.set_value('initial_conditions/cp/cpp_code', 'x[0] < 0.5 ? 0.0 : 1.0')
        sim.test_surf_normal = [1, 0, 0]
        sim.test_coord_index = 0
    elif surf_normal == 'y':
        sim.input.set_value('initial_conditions/cp/cpp_code', 'x[1] < 0.5 ? 1.0 : 0.0')
        sim.test_surf_normal = [0, -1, 0]
        sim.test_coord_index = 1
    elif surf_normal == 'z':
        sim.input.set_value('initial_conditions/cp/cpp_code', 'x[2] < 0.5 ? 1.0 : 0.0')
        sim.test_surf_normal = [0, 0, -1]
        sim.test_coord_index = 2

    sim.input.set_value('multiphase_solver/polynomial_degree_colour', vof_deg)

    # Run any custom code
    if modifier is not None:
        modifier(sim)

    sim.log.setup()
    sim.setup()
    sim.data['c'].assign(sim.data['cp'])
    return sim


@pytest.mark.parametrize('case_name', ['DG0_2D_y', 'DG0_3D_z', 'DG0_2D_x'])
def test_surface_locator(case_name):
    from ocellaris.probes.free_surface_locator import get_free_surface_locator

    counter = 0

    def hook():
        nonlocal counter
        counter += 1

    # Get a free surface locator
    vof_sim = mk_vof_sim(case_name)
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


@pytest.mark.parametrize(
    'case_name,extra_tall',
    [('DG0_2D_y', False), ('DG0_3D_z', False), ('DG0_2D_x', False), ('DG0_3D_z', True)],
)
def test_level_set_view(case_name, extra_tall):
    """
    Test creating LevelSetView functions in 2D and 3D with varying
    interface normals.

    The extra_tall parameter creates a narrow domain such that (when
    running on a size 3 or larger MPI domain) no surface cell distances
    will propagate to the rank containing the upper parts of the mesh.
    This is assuming that the mesh distribution algorithm does not do
    something stupid.
    """
    modifier = None

    coordmax = 1.0
    error_lim = 0.02
    if extra_tall:
        # Stretch the domain
        coordmax = 2.0
        error_lim = 0.05

        def modifier(sim):
            sim.input.set_value('mesh/endz', coordmax)
            sim.input.set_value('mesh/Nx', 2)
            sim.input.set_value('mesh/Ny', 2)
            sim.input.set_value('mesh/Nz', 20)

    # Get the level set view and force an update
    vof_sim = mk_vof_sim(case_name, modifier)
    lsv = vof_sim.multi_phase_model.get_level_set_view()
    vof_sim.hooks.run_custom_hook('MultiPhaseModelUpdated')

    # Computed and expected functions
    lsf = lsv.level_set_function
    V = lsf.function_space()
    cpp = 'std::abs(0.5 - x[%d])' % vof_sim.test_coord_index
    expected = dolfin.Expression(cpp, element=V.ufl_element())

    # Get the computed and expected value arrays
    arr = lsf.vector().get_local()
    arr_expected = dolfin.interpolate(expected, V).vector().get_local()
    print(arr)
    print(arr_expected)

    # Check that there are a significant number of values that are not the
    # default distance 1e100 away
    occurences_local = sum(1 if v < 1e100 else 0 for v in arr)
    occurences = dolfin.MPI.sum(dolfin.MPI.comm_world, float(occurences_local))
    assert occurences / V.dim() > 0.4  # Should be OK even with MPI.size == 2

    # Ignore 1e100 values when computing the error to enable checking the
    # error even when a rank has no free surface in its domain or in any
    # of its neighbour's domains such that all its values will be 1e10
    arr_expected = dolfin.interpolate(expected, V).vector().get_local()
    sieve = numpy.not_equal(arr, coordmax - 0.5)

    # The LSV is only a quasi distance function due to following the
    # edges of the mesh which exagerates distances a bit. Scale back to
    # max abs value to coordmax - 0.5 for comparing with the expected sol.
    amax = dolfin.MPI.max(dolfin.MPI.comm_world, float(abs(arr).max()))
    scale = (coordmax - 0.5) / amax

    diff = arr_expected[sieve] - scale * arr[sieve]
    error_local = diff.dot(diff)
    error = dolfin.MPI.sum(dolfin.MPI.comm_world, float(error_local))
    error = error ** 0.5 / occurences

    # 2D plot used when debugging
    if '2D' in vof_sim.test_name and False:
        from matplotlib import pyplot

        fig = pyplot.figure()
        c = dolfin.plot(lsf)
        pyplot.colorbar(c)
        pyplot.title(vof_sim.test_name + ' - error %g' % error)
        fig.savefig(vof_sim.test_name + '.png')

    # 3D plot used when debugging
    if extra_tall and True:
        with dolfin.XDMFFile(vof_sim.test_name + '.xdmf') as xdmf:
            xdmf.write(lsf)

    assert error < error_lim
