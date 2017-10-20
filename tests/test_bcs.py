import numpy
import dolfin
from ocellaris import Simulation, setup_simulation, run_simulation
from dolfin_utils.test import skip_in_parallel
from poisson_solver import BASE_INPUT
import pytest


@skip_in_parallel
@pytest.mark.parametrize("method", ['const', 'py_eval', 'py_exec', 'cpp'])
def test_dirichlet_bcs_scalar_constant_value(method):
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
    return (phi.replace('omega', omega),
            dphidx.replace('omega', omega),
            dphidy.replace('omega', omega),
            f.replace('omega', omega))


@skip_in_parallel
@pytest.mark.parametrize("method", ['py_eval', 'py_exec', 'cpp'])
def test_dirichlet_bcs_scalar_mms(method):
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
    relative_error = phidiff/analytical
    print('RELATIVE ERROR IS %.3f' % relative_error)
    assert relative_error < 0.074


@skip_in_parallel
@pytest.mark.parametrize("method", ['py_eval', 'py_exec', 'cpp'])
def test_neumann_bcs_scalar_mms(method):
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    sim.input.set_value('mesh/Nx', 20)
    sim.input.set_value('mesh/Ny', 20)
    sim.input.set_value('boundary_conditions', [{}, {}])
    sim.input.set_value('boundary_conditions/0/name', 'vertical walls')
    sim.input.set_value('boundary_conditions/0/selector', 'code')
    sim.input.set_value('boundary_conditions/0/inside_code',
                        'on_boundary and (x[0] < 1e-6 or x[0] > 1 - 1e-6)')
    sim.input.set_value('boundary_conditions/1/name', 'horizontal walls')
    sim.input.set_value('boundary_conditions/1/selector', 'code')
    sim.input.set_value('boundary_conditions/1/inside_code',
                        'on_boundary and (x[1] < 1e-6 or x[1] > 1 - 1e-6)')
    
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
    phih.vector()[:] += phia.vector()[0] - phih.vector()[0]
    phih.vector().apply('insert')
    
    # Compute relative error and check that it is reasonable
    phidiff = dolfin.errornorm(phi, phih)
    analytical = dolfin.norm(phia)
    relative_error = phidiff/analytical
    print('RELATIVE ERROR IS %.3f' % relative_error)
    assert relative_error < 0.055
