import dolfin
from ocellaris import Simulation, setup_simulation
from ocellaris.solver_parts.slope_limiter import SlopeLimiter
from poisson_solver import BASE_INPUT
from helpers import comm_self_to_comm_world
import pytest


def mk_limiter(degree, dim, use_cpp, comm_self=False):
    dolfin.parameters['ghost_mode'] = 'shared_vertex'

    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    if dim == 2:
        sim.input.set_value('mesh/type', 'Rectangle')
        sim.input.set_value('mesh/Nx', 10)
        sim.input.set_value('mesh/Ny', 10)
        cpp = 'A + A*sin(B*pi*x[0])*sin(B*pi*x[1])'
    else:
        sim.input.set_value('mesh/type', 'Box')
        sim.input.set_value('mesh/Nx', 5)
        sim.input.set_value('mesh/Ny', 5)
        sim.input.set_value('mesh/Nz', 5)
        cpp = 'A + A*sin(B*pi*x[0])*sin(B*pi*x[1])*sin(B*pi*x[2])'

    # Create simulation with mesh and the phi function
    sim.input.set_value('mesh/mpi_comm', 'SELF' if comm_self else 'WORLD')
    sim.input.set_value('solver/polynomial_degree', degree)
    sim.input.set_value('output/stdout_enabled', False)
    setup_simulation(sim)
    phi = sim.data['phi']
    V = phi.function_space()

    # Create a phi field with some jumps
    e = dolfin.Expression(cpp, element=V.ufl_element(), A=0.5, B=2.0)
    phi.interpolate(e)
    arr = phi.vector().get_local()
    arr[arr > 0.8] = 2 * 0.8 - arr[arr > 0.8]
    arr[arr < 0.2] = 0.5
    phi.vector().set_local(arr)
    phi.vector().apply('insert')

    # Create slope limiter
    sim.input.set_value('slope_limiter/phi/method', 'HierarchicalTaylor')
    sim.input.set_value('slope_limiter/phi/use_cpp', use_cpp)
    sim.input.set_value('slope_limiter/phi/skip_boundaries', [])
    lim = SlopeLimiter(sim, 'phi', phi)

    if True:
        from matplotlib import pyplot
        pyplot.figure()
        patches = dolfin.plot(phi)
        pyplot.colorbar(patches)
        pyplot.savefig('ht_lim_phi.png')

    return phi, lim


@pytest.mark.parametrize("degree,dim,test_mpi", [(1, 2, True), (2, 3, True),
                                                 (1, 2, False), (1, 3, False),
                                                 (2, 2, False), (2, 3, False)][:2])
def test_htlim_cpp_vs_py(degree, dim, test_mpi):
    """
    Apply the Hierarchical Taylor slope limiter with and without C++
    to verify the C++ implementation correctness.

    FIXME: Python/C++ comparisons disabled ([:2] above) due to missing BC impl. in Python

    Also run parallel vs serial and check equivalence
    """
    # Run HT lim in Python and C++
    phis = []
    for use_cpp in (False, True):
        comm_self = False
        if test_mpi:
            # Run Python code in serial, C++ code in parallel
            comm_self = not use_cpp
            use_cpp = True  # FIXME: HACK due to to missing BC impl. in Python

        # Get unlimited phi
        phi, lim = mk_limiter(degree, dim, use_cpp, comm_self)
        n1 = phi.vector().norm('l2')

        # Limit phi
        lim.run()
        n2 = phi.vector().norm('l2')

        # Sanity checks
        print(repr(lim), use_cpp, n1, n2)
        assert 'HierarchicalTaylor' in repr(lim)
        assert n1 != n2
        assert not comm_self or phi.function_space().mesh().mpi_comm().size == 1
        phis.append(phi)

    # Check results
    p0 = phis[0].vector()
    p1 = phis[1].vector()
    if test_mpi:
        p0 = comm_self_to_comm_world(phis[0], phis[1]).vector()

    diff = p0.copy()
    diff.axpy(-1, p1)
    p0n = p0.norm('l2')
    dn = diff.norm('l2')
    print(degree, dim, test_mpi, p0n, dn)
    assert 50 > p0n > 5  # not close to zero or very large
    assert dn < p0n / 1e15  # relative diff
