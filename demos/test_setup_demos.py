import os
from ocellaris import Simulation, setup_simulation
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
    skip = {'linear_sloshing.inp', 'falling_sphere.inp', 'dead_water_2D.inp'}
    if demo_inp_file in skip:
        return pytest.skip()
    
    monkeypatch.chdir(DEMODIR)
    sim = Simulation()
    sim.input.read_yaml(demo_inp_file)
    success = setup_simulation(sim)
    assert success
