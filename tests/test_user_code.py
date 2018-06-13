import sys, os
from ocellaris import Simulation, setup_simulation


TEST_DIR = os.path.dirname(__file__)


# Dummy values to make setup_simulation() run without errors
DUMMY = """
ocellaris:
    type: input
    version: 1.0
solver: {type: AnalyticalSolution}
mesh: {type: Rectangle, Nx: 4, Ny: 4}
boundary_conditions: []
time: {dt: 1.0}
physical_properties: {nu: 1.0}
output: {log_enabled: no}
"""


INPUT_USER_CONSTANTS = (
    """
user_code:
    constants:
        A: 21

ref: py$ A*2.0
"""
    + DUMMY
)


def test_user_constants():
    sim = Simulation()
    sim.input.read_yaml(yaml_string=INPUT_USER_CONSTANTS)
    success = setup_simulation(sim)
    assert success

    assert sim.input.get_value('ref') == 42.0


INPUT_MODULE_IMPORT = (
    """
user_code:
    python_path:
    -   %s/scripts
    modules:
    -   plot_reports
"""
    + DUMMY
)


def test_import_module():
    # A randomly selected script that does nothing at import time
    dummy_mod = 'plot_reports'
    assert dummy_mod not in sys.modules

    # Account for current working directory
    inp = INPUT_MODULE_IMPORT % os.path.join(TEST_DIR, '..')

    sim = Simulation()
    sim.input.read_yaml(yaml_string=inp)
    success = setup_simulation(sim)
    assert success

    assert dummy_mod in sys.modules
    sys.modules.pop(dummy_mod)


INPUT_USER_CODE = (
    """
user_code:
    constants:
        Q: 1.0
    code: |
        import sys
        sys.modules['__DUMMY__'] = 1
        assert simulation is not None
        assert Q == 1.0
"""
    + DUMMY
)


def test_use_code():
    dummy_mod = '__DUMMY__'
    assert dummy_mod not in sys.modules

    sim = Simulation()
    sim.input.read_yaml(yaml_string=INPUT_USER_CODE)
    success = setup_simulation(sim)
    assert success

    assert dummy_mod in sys.modules
    sys.modules.pop(dummy_mod)
