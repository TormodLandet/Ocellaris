import os, sys, subprocess
import pytest


DEMODIR = os.path.abspath(os.path.dirname(__file__))
def pytest_generate_tests(metafunc):
    """
    Setup input files to test (all *.inp files in the same
    directory as this file
    """
    inpfiles = []
    for fn in os.listdir(DEMODIR):
        if fn.endswith('.inp'):
            inpfiles.append(fn)
    metafunc.parametrize("demo_inp_file", inpfiles)


def test_setup_demo(demo_inp_file, monkeypatch):
    "Run setup on demo input file"
    skip = {'linear_sloshing.inp', 'falling_sphere.inp'}
    if demo_inp_file in skip:
        return pytest.skip()
    
    slow = {'internal_soliton.inp', 'dead_water_2D.inp', 'dead_water_3D.inp'}
    if demo_inp_file in slow and os.environ.get('OCELLARIS_RUN_SLOW_TEST') != '1':
        raise pytest.skip('Skipping slow test')
    
    monkeypatch.chdir(DEMODIR)
    
    # run python3 -m ocellaris DEMOFILE --set-inp time/tmax=0
    interpreter = sys.executable
    mainfile = os.path.join(DEMODIR, '..', 'ocellaris', '__main__.py')
    cmd = [interpreter, mainfile, demo_inp_file, '--set-inp', 'time/tmax=0']
    subprocess.check_call(cmd)
    
    # Run as script (can lead to inter-demo disturbances with dolfin globals etc)
    #from ocellaris import Simulation, setup_simulation
    #sim = Simulation()
    #sim.input.read_yaml(demo_inp_file)
    #success = setup_simulation(sim)
    #assert success
