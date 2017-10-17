import numpy
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
    assert numpy.linalg.norm(p - 1.0) < 1e-12
