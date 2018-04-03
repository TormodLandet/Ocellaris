from helpers import get_test_file_name
from ocellaris import Simulation


def test_base_input():
    fn = get_test_file_name('base.inp')
    sim = Simulation()
    sim.input.read_yaml(fn)
    
    gv = lambda k, t: sim.input.get_value(k, required_type=t) 
    assert gv('metadata/author', 'string') == 'Tormod Landet'
    assert gv('metadata/description', 'string') == 'NoDescription'
    
    assert tuple(gv('some_vals/bools', 'list(int)')) == (1, 1, 0, 0)
    assert tuple(gv('some_vals/floats', 'list(float)')) == (1.1, 2, 3.0e3)
    assert gv('some_vals/computed', 'float') == 2.0


def test_child_input():
    fn = get_test_file_name('child.inp')
    sim = Simulation()
    sim.input.read_yaml(fn)
    
    gv = lambda k, t: sim.input.get_value(k, required_type=t) 
    assert gv('metadata/author', 'string') == 'Tormod Landet'
    #assert gv('metadata/date', 'date') == '2018-04-03'
    assert gv('metadata/description', 'string') == 'ThisIsTheDescription'
    
    assert gv('some_vals/computed', 'float') == 2.0
    assert gv('test_py_val/C', 'float') == 8.0
    