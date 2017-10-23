"""
Test restart file handling and general IO
"""
import numpy
from ocellaris import Simulation, setup_simulation
from dolfin_utils.test import skip_in_parallel
from poisson_solver import BASE_INPUT


@skip_in_parallel
def test_restart_file_io(tmpdir):
    prefix = str(tmpdir.mkdir("test_restart_file_io").join('ocellaris'))
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    sim.input.set_value('output/prefix', prefix)
    sim.input.set_value('time/tstart', 42.0)
    setup_simulation(sim)
    
    # Fill in the phi function
    phi = sim.data['phi']
    phi_arr = phi.vector().get_local()
    phi_arr[:] = numpy.random.rand(*phi_arr.shape)
    phi.vector().set_local(phi_arr)
    phi.vector().apply('insert')
    
    # Save restart file
    filename = sim.io.write_restart_file()
    assert filename.startswith(prefix)
    
    # Load input from restart file
    sim2 = Simulation()
    sim2.io.load_restart_file_input(filename)
    assert sim2.input.get_value('time/tstart') == 42.0
    assert str(sim.input) == str(sim2.input)
    
    # Load phi from restart file
    setup_simulation(sim2)
    sim2.io.load_restart_file_results(filename)
    phi2 = sim2.data['phi']
    phi2_arr = phi2.vector().get_local()
    assert all(phi_arr == phi2_arr)
    
    assert sim.data['mesh'].hash() == sim2.data['mesh'].hash() 
