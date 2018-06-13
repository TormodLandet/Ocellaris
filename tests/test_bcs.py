"""
Test solving the Poisson equation with different boundary conditions
This also test the whole Ocellaris solver framework with input file
reading and the setup/run functionality
"""
import numpy
import dolfin
from ocellaris import Simulation, setup_simulation, run_simulation
from poisson_solver import BASE_INPUT
import pytest


@pytest.mark.parametrize("method", ['const', 'py_eval', 'py_exec', 'cpp'])
def test_dirichlet_bcs_scalar_constant_value(method):
    "Test inhomogenous Dirichlet BCs using a Poisson solver"
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    sim.input.set_value('boundary_conditions', [{}])
    sim.input.set_value('boundary_conditions/0/name', 'all walls')
    sim.input.set_value('boundary_conditions/0/selector', 'code')
    sim.input.set_value('boundary_conditions/0/inside_code', 'on_boundary')

    if method == 'const':
        sim.input.set_value('boundary_conditions/0/phi/type', 'ConstantValue')
        sim.input.set_value('boundary_conditions/0/phi/value', 1.0)
    elif method == 'py_eval':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CodedValue')
        sim.input.set_value('boundary_conditions/0/phi/code', '1.0')
    elif method == 'py_exec':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CodedValue')
        sim.input.set_value('boundary_conditions/0/phi/code', 'value[0] = 1.0')
    elif method == 'cpp':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CppCodedValue')
        sim.input.set_value('boundary_conditions/0/phi/cpp_code', '1.0')

    setup_simulation(sim)
    run_simulation(sim)

    p = sim.data['phi'].vector().get_local()
    assert numpy.linalg.norm(p - 1.0) < 1e-8


def mms_case(omega='1.23'):
    """
    A manufactured solution of

      -∇⋅∇φ = f
    
    Returns C++ code for expressions
    """
    phi = 'sin(omega*pi*x[0])*sin(omega*pi*x[1])'
    dphidx = 'pi*omega*cos(omega*pi*x[0])*sin(omega*pi*x[1])'
    dphidy = 'pi*omega*sin(omega*pi*x[0])*cos(omega*pi*x[1])'
    f = '2*pi*pi*omega*omega*sin(omega*pi*x[0])*sin(omega*pi*x[1])'
    return (
        phi.replace('omega', omega),
        dphidx.replace('omega', omega),
        dphidy.replace('omega', omega),
        f.replace('omega', omega),
    )


@pytest.mark.parametrize("method", ['py_eval', 'py_exec', 'cpp'])
def test_dirichlet_bcs_scalar_mms(method):
    "Test inhomogenous coded Dirichlet BCs using a Poisson solver"
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    sim.input.set_value('boundary_conditions', [{}])
    sim.input.set_value('boundary_conditions/0/name', 'all walls')
    sim.input.set_value('boundary_conditions/0/selector', 'code')
    sim.input.set_value('boundary_conditions/0/inside_code', 'on_boundary')

    # Get analytical expressions
    cphi, _, _, cf = mms_case()
    sim.input.set_value('solver/source', cf)

    # Setup the boundary conditions to test
    if method == 'py_eval':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CodedValue')
        sim.input.set_value('boundary_conditions/0/phi/code', cphi)
    elif method == 'py_exec':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CodedValue')
        sim.input.set_value('boundary_conditions/0/phi/code', 'value[0] = ' + cphi)
    elif method == 'cpp':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CppCodedValue')
        sim.input.set_value('boundary_conditions/0/phi/cpp_code', cphi)

    # Run Ocellaris
    setup_simulation(sim)
    run_simulation(sim)

    # The numeric (phih) and analytic (phia) solution functions
    Vphi = sim.data['Vphi']
    phi = dolfin.Expression(cphi, degree=5)
    phih = sim.data['phi']
    phia = dolfin.interpolate(phi, Vphi)

    # Compute relative error and check that it is reasonable
    phidiff = dolfin.errornorm(phi, phih)
    analytical = dolfin.norm(phia)
    relative_error = phidiff / analytical
    print('RELATIVE ERROR IS %.3f' % relative_error)
    assert relative_error < 0.074


def correct_constant_offset(sim, phih, phia):
    """
    Correct the constant offset due to how the null space is handled
    """
    if sim.rank == 0:
        mod = phia.vector()[0] - phih.vector()[0]
        mod = dolfin.MPI.max(dolfin.MPI.comm_world, mod)
    else:
        mod = dolfin.MPI.max(dolfin.MPI.comm_world, -1.0e100)
    phih.vector()[:] += mod
    phih.vector().apply('insert')


@pytest.mark.parametrize("method", ['py_eval', 'py_exec', 'cpp'])
def test_neumann_bcs_scalar_mms(method):
    "Test pure Neumann BCs using a Poisson solver"
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    sim.input.set_value('mesh/Nx', 20)
    sim.input.set_value('mesh/Ny', 20)
    sim.input.set_value('boundary_conditions', [{}, {}])
    sim.input.set_value('boundary_conditions/0/name', 'vertical walls')
    sim.input.set_value('boundary_conditions/0/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/0/inside_code',
        'on_boundary and (x[0] < 1e-6 or x[0] > 1 - 1e-6)',
    )
    sim.input.set_value('boundary_conditions/1/name', 'horizontal walls')
    sim.input.set_value('boundary_conditions/1/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/1/inside_code',
        'on_boundary and (x[1] < 1e-6 or x[1] > 1 - 1e-6)',
    )

    # Get analytical expressions
    cphi, cphix, cphiy, cf = mms_case()
    sim.input.set_value('solver/source', cf)

    if method.startswith('py'):
        cphix = '(-1.0 if x[0] < 0.5 else 1.0) * ' + cphix
        cphiy = '(-1.0 if x[1] < 0.5 else 1.0) * ' + cphiy
    else:
        cphix = '(x[0] < 0.5 ? -1.0 : 1.0) * ' + cphix
        cphiy = '(x[1] < 0.5 ? -1.0 : 1.0) * ' + cphiy

    # Setup the boundary conditions to test
    if method == 'py_eval':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CodedGradient')
        sim.input.set_value('boundary_conditions/0/phi/code', cphix)
        sim.input.set_value('boundary_conditions/1/phi/type', 'CodedGradient')
        sim.input.set_value('boundary_conditions/1/phi/code', cphiy)
    elif method == 'py_exec':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CodedGradient')
        sim.input.set_value('boundary_conditions/0/phi/code', 'value[0] = ' + cphix)
        sim.input.set_value('boundary_conditions/1/phi/type', 'CodedGradient')
        sim.input.set_value('boundary_conditions/1/phi/code', 'value[0] = ' + cphiy)
    elif method == 'cpp':
        sim.input.set_value('boundary_conditions/0/phi/type', 'CppCodedGradient')
        sim.input.set_value('boundary_conditions/0/phi/cpp_code', cphix)
        sim.input.set_value('boundary_conditions/1/phi/type', 'CppCodedGradient')
        sim.input.set_value('boundary_conditions/1/phi/cpp_code', cphiy)

    # Run Ocellaris
    setup_simulation(sim)
    run_simulation(sim)

    # The numeric (phih) and analytic (phia) solution functions
    Vphi = sim.data['Vphi']
    phi = dolfin.Expression(cphi, degree=5)
    phih = sim.data['phi']
    phia = dolfin.interpolate(phi, Vphi)

    # Correct the constant offset due to how the null space is handled
    correct_constant_offset(sim, phih, phia)

    # Compute relative error and check that it is reasonable
    phidiff = dolfin.errornorm(phi, phih)
    analytical = dolfin.norm(phia)
    relative_error = phidiff / analytical
    print('RELATIVE ERROR IS %.3f' % relative_error)
    assert relative_error < 0.055


@pytest.mark.parametrize(
    "bcs,b", [('robin', 1), ('robin', 0.1), ('robin', 0.01), ('neumann', None)]
)
def test_robin_bcs_scalar_mms(bcs, b):
    """
    Test Robin BCs using a Poisson solver to solve
    
      -∇⋅∇φ = f
    
    where φ = 1 + x and hence f = 0. We use Neumann
    BCs n⋅∇φ = 0 on the horizontal walls and Robin
    BCs on the vertical walls
    """
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    # Create boundary regions
    sim.input.set_value('boundary_conditions', [{}, {}, {}])
    sim.input.set_value('boundary_conditions/0/name', 'vertical wall x=0')
    sim.input.set_value('boundary_conditions/0/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/0/inside_code',
        'on_boundary and (x[0] < 1e-6 or x[0] > 1 - 1e-6)',
    )
    sim.input.set_value('boundary_conditions/1/name', 'vertical walls x=1')
    sim.input.set_value('boundary_conditions/1/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/1/inside_code', 'on_boundary and x[0] > 1 - 1e-6'
    )
    sim.input.set_value('boundary_conditions/2/name', 'horizontal walls')
    sim.input.set_value('boundary_conditions/2/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/2/inside_code',
        'on_boundary and (x[1] < 1e-6 or x[1] > 1 - 1e-6)',
    )

    # Setup the boundary conditions to test
    if bcs == 'robin':
        sim.input.set_value('boundary_conditions/0/phi/type', 'ConstantRobin')
        sim.input.set_value('boundary_conditions/0/phi/blend', b)
        sim.input.set_value('boundary_conditions/0/phi/dval', 1.0)
        sim.input.set_value('boundary_conditions/0/phi/nval', -1.0)
        sim.input.set_value('boundary_conditions/1/phi/type', 'ConstantRobin')
        sim.input.set_value('boundary_conditions/1/phi/blend', b)
        sim.input.set_value('boundary_conditions/1/phi/dval', 2.0)
        sim.input.set_value('boundary_conditions/1/phi/nval', 1.0)
        sim.input.set_value('boundary_conditions/2/phi/type', 'ConstantGradient')
        sim.input.set_value('boundary_conditions/2/phi/value', 0.0)
    elif bcs == 'neumann':
        sim.input.set_value('boundary_conditions/0/phi/type', 'ConstantGradient')
        sim.input.set_value('boundary_conditions/0/phi/value', -1.0)
        sim.input.set_value('boundary_conditions/1/phi/type', 'ConstantGradient')
        sim.input.set_value('boundary_conditions/1/phi/value', 1.0)

    # RHS
    sim.input.set_value(
        'solver/source', '2*pi*pi*(pow(sin(x[0]*pi), 2) - pow(cos(x[0]*pi), 2))'
    )

    # Run Ocellaris
    setup_simulation(sim)
    run_simulation(sim)

    # The numeric (phih) and analytic (phia) solution functions
    cphi = '1.0 + x[0] + pow(sin(x[0]*pi), 2)'
    Vphi = sim.data['Vphi']
    phi = dolfin.Expression(cphi, degree=5)
    phih = sim.data['phi']
    phia = dolfin.interpolate(phi, Vphi)

    # Correct the constant offset due to how the null space is handledFalse
    if bcs == 'neumann':
        correct_constant_offset(sim, phih, phia)

    # Plot to file for debugging
    # debug_phi_plot(phia, phih, 'test_robin_bcs_scalar_mms_%s.png' % bcs)

    # Compute relative error and check that it is reasonable
    phidiff = dolfin.errornorm(phi, phih)
    analytical = dolfin.norm(phia)
    relative_error = phidiff / analytical
    print('RELATIVE ERROR IS %.3f' % relative_error)
    assert relative_error < 0.015  # Expect 0.0139 with Robin


@pytest.mark.parametrize(
    "method,slip_length",
    [
        ('Constant', 1000.0),
        ('Constant', 1.0),
        ('Constant', 0.001),
        ('Constant', 0.0),
        ('C++', 1.0),
        ('Interface', 1.0),
    ],
)
def test_slip_length_robin_bcs_scalar_mms(slip_length, method):
    """
    Test slip length Robin BCs using a Poisson solver to solve
    
      -∇⋅∇φ = f
    
    where φ = (-6x² + 6x + 6𝛿)/(6𝛿  - 1) and hence f = 12/(6𝛿  - 1). 
    
    We use Neumann BCs n⋅∇φ = 0 on the horizontal walls and Navier's
    slip length boundary condition on the vertical walls. The selectedFalse
    analytical solution is such that for any slip length 𝛿 the average
    value of φ is 1.0
    
    This mimics a flow profile going vertically in a 1.0 wide channel
    """
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    # Create boundary regions
    sim.input.set_value('boundary_conditions', [{}, {}])
    sim.input.set_value('boundary_conditions/0/name', 'vertical walls')
    sim.input.set_value('boundary_conditions/0/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/0/inside_code',
        'on_boundary and (x[0] < 1e-6 or x[0] > 1 - 1e-6)',
    )
    sim.input.set_value('boundary_conditions/1/name', 'horizontal walls')
    sim.input.set_value('boundary_conditions/1/selector', 'code')
    sim.input.set_value(
        'boundary_conditions/1/inside_code',
        'on_boundary and (x[1] < 1e-6 or x[1] > 1 - 1e-6)',
    )

    # Vertical wall BCs
    if method == 'Constant':
        sim.input.set_value('boundary_conditions/0/phi/type', 'SlipLength')
        sim.input.set_value('boundary_conditions/0/phi/slip_length', slip_length)
    elif method == 'C++':
        sim.input.set_value('boundary_conditions/0/phi/type', 'SlipLength')
        sim.input.set_value('boundary_conditions/0/phi/slip_length', repr(slip_length))
    else:
        sim.input.set_value('boundary_conditions/0/phi/type', 'InterfaceSlipLength')
        sim.input.set_value('boundary_conditions/0/phi/slip_length', slip_length)
        sim.input.set_value('boundary_conditions/0/phi/slip_factor_distance', 0.1)

        # Create dummy scalar field
        sim.input.set_value('multiphase_solver/type', 'BlendedAlgebraicVOF')
        sim.input.set_value('multiphase_solver/function_space_colour', 'DG')
        sim.input.set_value('multiphase_solver/polynomial_degree_colour', 0)
        sim.input.set_value(
            'initial_conditions/cp/cpp_code', '0.5'
        )  # The interface is everywhere

        sim.input.set_value('physical_properties/rho0', 1.0)
        sim.input.set_value('physical_properties/rho1', 1.0)
        sim.input.set_value('physical_properties/nu0', 1.0)
        sim.input.set_value('physical_properties/nu1', 1.0)

        sim.input.set_value('multiphase_solver/project_uconv_dgt0', False)
        sim.data['u_conv'] = dolfin.as_vector([0, 0])
        sim.data['dt'] = dolfin.Constant(1.0)

    # Horizontal wall BCs
    sim.input.set_value('boundary_conditions/1/phi/type', 'ConstantGradient')
    sim.input.set_value('boundary_conditions/1/phi/value', 0.0)

    # RHS
    sim.input.set_value(
        'solver/source', '12/(6*𝛿  - 1.0)'.replace('𝛿', repr(slip_length))
    )

    # Run Ocellaris
    setup_simulation(sim)
    run_simulation(sim)

    # The numeric (phih) and analytic (phia) solution functions
    cphi = '(-6*x[0]*x[0] + 6*x[0] + 6*𝛿)/(6*𝛿  - 1.0)'.replace('𝛿', repr(slip_length))
    Vphi = sim.data['Vphi']
    phi = dolfin.Expression(cphi, degree=5)
    phih = sim.data['phi']
    phia = dolfin.interpolate(phi, Vphi)

    # Plot to file for debugging
    # debug_phi_plot(phi, phia, phih, 'test_slip_length_bcs_scalar_mms_%g.png' % slip_length)

    # Compute relative error and check that it is reasonable
    phidiff = dolfin.errornorm(phi, phih)
    analytical = dolfin.norm(phia)
    relative_error = phidiff / analytical
    print('RELATIVE ERROR IS %.4f for 𝛿=%r' % (relative_error, slip_length))

    assert relative_error < 0.0099  # Expect 0.0097 for 𝛿=0.001


@pytest.mark.parametrize("N", [4, 6])
def test_slip_length_boundary_level_set_intersector(N):
    from ocellaris.solver_parts.boundary_conditions.slip_length import (
        BoundaryLevelSetIntersector
    )

    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)

    # Create simulation with a scalar field that has a level set at x[1] == 0.5
    sim.input.set_value('mesh/Nx', N)
    sim.input.set_value('mesh/Ny', N)
    sim.input.set_value('multiphase_solver/type', 'BlendedAlgebraicVOF')
    sim.input.set_value('multiphase_solver/function_space_colour', 'DG')
    sim.input.set_value('multiphase_solver/polynomial_degree_colour', 0)
    sim.input.set_value('initial_conditions/cp/cpp_code', 'x[1] < 0.5 ? 1.0 : 0.0')
    sim.input.set_value('physical_properties/rho0', 1.0)
    sim.input.set_value('physical_properties/rho1', 1.0)
    sim.input.set_value('physical_properties/nu0', 1.0)
    sim.input.set_value('physical_properties/nu1', 1.0)
    setup_simulation(sim)

    intersector = BoundaryLevelSetIntersector(sim, sim.data['cp'], 0.5)
    intersections = intersector.get()

    for pos in intersections:
        print('Intersection: %r' % pos)

    assert len(intersections) == 2
    for pos in intersections:
        assert numpy.allclose(pos, [0, 0.5]) or numpy.allclose(pos, [1, 0.5])


def debug_phi_plot(phi, phia, phih, plotname, compare_historical=False, **plotargs):
    """
    Only used while developing the tests, left here in case more tests are
    added that need some debug plotting as well
    """
    diff = phia.copy(deepcopy=True)
    diff.vector().axpy(-1, phih.vector())
    diff.vector().apply('insert')

    na = dolfin.norm(phia)
    nh = dolfin.norm(phih)
    nd = dolfin.norm(diff)
    nd2 = dolfin.errornorm(phi, phih)

    from matplotlib import pyplot

    pyplot.figure(figsize=(10, 20))

    pyplot.subplot(311)
    c = dolfin.plot(phia, **plotargs)
    pyplot.colorbar(c)
    pyplot.title('Analytical (norm: %g)' % na)

    pyplot.subplot(312)
    c = dolfin.plot(phih, **plotargs)
    pyplot.colorbar(c)
    pyplot.title('Numerical (norm: %g)' % nh)

    pyplot.subplot(313)
    c = dolfin.plot(diff, **plotargs)
    pyplot.colorbar(c)
    pyplot.title('Diff (norm: %g, errornorm: %g, rel: %g)' % (nd, nd2, nd2 / na))

    pyplot.tight_layout()
    pyplot.savefig(plotname)
    pyplot.close()

    if not compare_historical:
        return
    hist_file = plotname + '.npy'

    comm = phih.function_space().mesh().mpi_comm()
    mpi_size = dolfin.MPI.size(comm)
    if mpi_size != 1:
        return

    import os

    arr = phih.vector().get_local()
    if not os.path.isfile(hist_file):
        numpy.save(hist_file, arr)
        return

    hist = numpy.load(hist_file)
    print('hist', hist[:5])
    print('arr ', arr[:5])
    print('diff', (hist - arr)[:5])
    print('maxabs', abs(hist - arr).max())

    assert abs(hist - arr).max() < 1e-12
