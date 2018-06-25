from dolfin import (
    UnitSquareMesh,
    UnitCubeMesh,
    FunctionSpace,
    Function,
    Expression,
    as_vector,
    MPI,
    parameters,
)
from ocellaris import Simulation
from ocellaris.solver_parts.convection import (
    get_convection_scheme,
    VelocityDGT0Projector,
)
from helpers import comm_self_to_comm_world
import pytest


def mk_scheme(N, Vname, Vorder, cpp_expr, expr_args, convection_inp, dim=2, comm=None):
    if comm is None:
        comm = MPI.comm_world

    parameters['ghost_mode'] = 'shared_vertex'
    if dim == 2:
        mesh = UnitSquareMesh(comm, N, N)
    else:
        mesh = UnitCubeMesh(comm, N, N, N)

    V = FunctionSpace(mesh, Vname, Vorder)
    C = Function(V)
    e = Expression(cpp_expr, element=V.ufl_element(), **expr_args)
    C.interpolate(e)

    D = Function(V)
    D.assign(C)

    sim = Simulation()
    sim.set_mesh(mesh)
    sim.data['constrained_domain'] = None
    sim.data['C'] = C
    for key, value in convection_inp.items():
        sim.input.set_value('convection/C/%s' % key, value)

    scheme_name = convection_inp['convection_scheme']
    return get_convection_scheme(scheme_name)(sim, 'C')


def mk_vel(sim, Vname, Vorder, cpp_exprs):
    mesh = sim.data['mesh']
    V = FunctionSpace(mesh, Vname, Vorder)

    vel = []
    for cpp in cpp_exprs:
        u = Function(V)
        u.interpolate(Expression(cpp, element=V.ufl_element()))
        vel.append(u)

    return as_vector(vel)


def mk_blending_factor(convection_inp, dim=2, comm=None):
    """
    Create convected function and convection scheme
    """
    if dim == 2:
        N = 10
        cpp = 'A + A*sin(B*pi*x[0])*sin(B*pi*x[1])'
    else:
        N = 5
        cpp = 'A + A*sin(B*pi*x[0])*sin(B*pi*x[1])*sin(B*pi*x[2])'

    scheme = mk_scheme(
        N,
        Vname='DG',
        Vorder=0,
        cpp_expr=cpp,
        expr_args={'A': 0.5, 'B': 2},
        convection_inp=convection_inp,
        dim=dim,
        comm=comm,
    )
    sim = scheme.simulation
    C = sim.data['C']
    beta = scheme.blending_function

    # Filter C to make it more like a free surface simulation
    Carr = C.vector().get_local()
    Carr[Carr > 0.8] = 1
    Carr[Carr < 0.2] = 0
    C.vector().set_local(Carr)
    C.vector().apply('insert')

    # Make convecting velocity
    vel = mk_vel(sim, 'DG', 2, ['1'] * dim)
    proj = VelocityDGT0Projector(sim, vel)
    proj.update()
    vel_dgt0 = proj.velocity

    # Run the convection scheme
    h = 1 / N
    dt = h * 0.3 / 2 ** 0.5
    if scheme.need_alpha_gradient:
        scheme.initialize_gradient()
        scheme.gradient_reconstructor.reconstruct()
    scheme.update(dt, vel_dgt0)

    if False:
        import dolfin
        from matplotlib import pyplot

        pyplot.figure()
        patches = dolfin.plot(C)
        pyplot.colorbar(patches)
        pyplot.savefig('test_conv_C.png')

        from ocellaris.utils.plotting_trace import plot_matplotlib_dgt

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        patches = plot_matplotlib_dgt(beta)
        fig.colorbar(patches, ax=ax)
        pyplot.savefig('test_conv_beta.png')

    return beta, sim


def test_hric_cpp():
    """
    Run HRIC with and without C++ implementation to  verify the C++ 
    implementation correctness
    """
    # regression test value
    ANS1 = 7.85326907878949600
    ANS2 = 0.42426406871193789

    cinp = {'convection_scheme': 'HRIC', 'HRIC_version': 'HRIC'}

    betas = []
    Cofs = []
    for use_cpp in (False, True):
        cinp['use_cpp'] = use_cpp
        beta, sim = mk_blending_factor(cinp)
        betas.append(beta)
        bn = beta.vector().norm('l2')

        Cof_max = sim.reporting.timestep_xy_reports['Cof_max'][-1]
        Cofs.append(Cof_max)

        eps1, eps2 = (1e-15, 1e-15) if MPI.comm_world.size == 1 else (0.02, 1e-14)
        print(use_cpp, repr(bn), bn - ANS1, eps1, repr(Cof_max), Cof_max - ANS2, eps2)
        assert abs(bn - ANS1) < eps1
        assert abs(Cof_max - ANS2) < eps2

    diff = betas[0].copy(deepcopy=True)
    diff.vector().axpy(-1, betas[1].vector())

    assert diff.vector().norm('l2') == 0


@pytest.mark.parametrize("dim,test_mpi", [(2, False), (3, False), (3, True)])
def test_hric_cpp_vs_py(dim, test_mpi):
    """
    Calculate the HRIC blending factor with and without C++ to verify the C++
    implementation correctness.
    
    Also run parallel vs serial and check equivalence
    """
    cinp = {'convection_scheme': 'HRIC', 'HRIC_version': 'HRIC'}

    # Produce gradients in Python and C++
    comm = None
    betas = []
    Cofs = []
    for use_cpp in (False, True):
        if test_mpi:
            # Run Python code in serial, C++ code in parallel
            comm = MPI.comm_world if use_cpp else MPI.comm_self
        cinp['use_cpp'] = use_cpp
        beta, sim = mk_blending_factor(cinp, dim, comm)
        Cof_max = sim.reporting.timestep_xy_reports['Cof_max'][-1]
        betas.append(beta)
        Cofs.append(Cof_max)

    # Check results
    assert 1 > Cofs[0] > 0.3
    assert abs(Cofs[0] - Cofs[1]) < 1e-12

    b0 = betas[0].vector()
    b1 = betas[1].vector()
    if test_mpi:
        b0 = comm_self_to_comm_world(betas[0], betas[1]).vector()

    diff = b0.copy()
    diff.axpy(-1, b1)
    b0n = b0.norm('l2')
    assert 15 > b0n > 5  # not close to zero or very large
    assert diff.norm('l2') < b0n / 1e15  # relative diff


@pytest.mark.parametrize("dim,test_mpi", [(2, False), (3, False), (3, True)])
def test_gradient_cpp_vs_py(dim, test_mpi):
    """
    Run the gradient reconstruction with and without C++ to verify the C++
    implementation correctness.
    
    Also run parallel vs serial and check equivalence
    """
    if dim == 2:
        N = 10
        cpp = 'A + A*sin(B*pi*x[0])*sin(B*pi*x[1])'
    else:
        N = 5
        cpp = 'A + A*sin(B*pi*x[0])*sin(B*pi*x[1])*sin(B*pi*x[2])'

    # Produce gradients in Python and C++
    comm = None
    schemes = []
    for use_cpp in (False, True):
        if test_mpi:
            # Run Python code in serial, C++ code in parallel
            comm = MPI.comm_world if use_cpp else MPI.comm_self

        convection_inp = {'use_cpp_gradient': use_cpp, 'convection_scheme': 'HRIC'}
        scheme = mk_scheme(
            N,
            Vname='DG',
            Vorder=0,
            cpp_expr=cpp,
            expr_args={'A': 0.5, 'B': 2},
            convection_inp=convection_inp,
            dim=dim,
            comm=comm,
        )

        scheme.initialize_gradient()
        scheme.gradient_reconstructor.reconstruct()
        schemes.append(scheme)

    # Check results
    for d in range(dim):
        g0 = schemes[0].gradient_reconstructor.gradient[d].vector()
        g1 = schemes[1].gradient_reconstructor.gradient[d].vector()

        if test_mpi:
            g0 = comm_self_to_comm_world(
                schemes[0].gradient_reconstructor.gradient[d],
                schemes[1].gradient_reconstructor.gradient[d],
            )
            g0 = g0.vector()

        diff = g0.copy()
        diff.axpy(-1, g1)
        g0n = g0.norm('l2')
        assert 20 > g0n > 10
        assert diff.norm('l2') < g0n / 1e15
