import dolfin
from ocellaris import Simulation
from ocellaris.solver_parts.fields.sharp_field import SharpField
from utils import check_vector_value_histogram, get_vector_value_histogram

INP = """
ocellaris:
    type: input
    version: 1.0

user_code:
    constants:
        A: 1

fields:
-   name: dens0
    type: SharpField
    variable_name: rho
    z: 0.5
    value_above: 1
    value_below: 1000
"""


def test_sharp_field_dg0():
    sim = Simulation()
    sim.input.read_yaml(yaml_string=INP)
    mesh = dolfin.UnitCubeMesh(2, 2, 2)
    sim.set_mesh(mesh)

    # Sharp jump at the cell boundaries
    field_inp = sim.input.get_value('fields/0', required_type='Input')
    f = SharpField(sim, field_inp).get_variable('rho')
    check_vector_value_histogram(f.vector(), {1: 24, 1000: 24}, round_digits=6)

    # Sharp jump in the middle of the cells
    field_inp['z'] = 0.25
    f = SharpField(sim, field_inp).get_variable('rho')
    hist = get_vector_value_histogram(f.vector())
    assert len(hist) == 4
    for k, v in hist.items():
        if round(k, 6) != 1:
            assert v == 8
        else:
            assert v == 24

    # Non-projected sharp jump between cells
    field_inp['local_projection'] = False
    field_inp['z'] = 0.50
    f = SharpField(sim, field_inp).get_variable('rho')
    check_vector_value_histogram(f.vector(), {1: 24, 1000: 24}, round_digits=6)

    # Non-projected sharp jump the middle of the cells
    field_inp['z'] = 0.25
    f = SharpField(sim, field_inp).get_variable('rho')
    hist = get_vector_value_histogram(f.vector())
    check_vector_value_histogram(f.vector(), {1: 40, 1000: 8}, round_digits=6)
