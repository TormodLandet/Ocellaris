import dolfin
from ocellaris import Simulation, setup_simulation
from ocellaris.solver_parts.boundary_conditions import \
    get_dof_region_marks, mark_cell_layers


BASE_INPUT = """
ocellaris:
    type: input
    version: 1.0

mesh:
    type: Rectangle
    Nx: 10
    Ny: 10

boundary_conditions:
-   name: not left
    selector: code
    inside_code: on_boundary
    p:
        type: ConstantValue
        value: 0.0
-   name: left
    selector: code
    inside_code: on_boundary and x[0] < 1e-6
    p:
        type: ConstantValue
        value: 1.0

solver:
    type: IPCS-A

output:
    log_enabled: no
    solution_properties: off
    xdmf_write_interval: 0
    save_restart_file_at_end: off

time:
    dt: 1.0
    tmax: 1.0
physical_properties:
    nu: 1.0
"""


def test_get_dof_region_marks():
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    setup_simulation(sim)
    Vp = sim.data['Vp']
    dofs_x = Vp.tabulate_dof_coordinates().reshape((-1, 2))

    drm = get_dof_region_marks(sim, Vp)
    num_in_region = [0, 0]
    for dof, marks in drm.items():
        x, y = dofs_x[dof]
        if x == 0:
            print(x, y, dof, marks)
        assert x == 0 or x == 1 or y == 0 or y == 1
        assert (0 in marks) == (x == 1 or y == 0 or y == 1)
        assert (1 in marks) == (x == 0)
        for mark in marks:
            num_in_region[mark] += 1

    assert mpi_int_sum(num_in_region[0]) == 30 + 29 * 2
    assert mpi_int_sum(num_in_region[1]) == 30


def test_mark_cell_layers():
    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT)
    setup_simulation(sim)
    mesh = sim.data['mesh']
    Vp = sim.data['Vp']

    # No named == all
    cells = mark_cell_layers(sim, Vp)
    for cid in cells:
        cell = dolfin.Cell(mesh, cid)
        mp = cell.midpoint()[:]
        assert mp[0] < 0.1 or mp[0] > 0.9 or mp[1] < 0.1 or mp[1] > 0.9
    assert mpi_int_sum(len(cells)) == 20 * 2 + 16 * 2

    # Test all
    cells = mark_cell_layers(sim, Vp, named_boundaries=['all'])
    for cid in cells:
        cell = dolfin.Cell(mesh, cid)
        mp = cell.midpoint()[:]
        assert mp[0] < 0.1 or mp[0] > 0.9 or mp[1] < 0.1 or mp[1] > 0.9
    assert mpi_int_sum(len(cells)) * 2 + 16 * 2

    # Test only left side
    cells = mark_cell_layers(sim, Vp, named_boundaries=['left'])
    for cid in cells:
        cell = dolfin.Cell(mesh, cid)
        mp = cell.midpoint()[:]
        assert mp[0] < 0.1
    assert mpi_int_sum(len(cells)) == 20

    cells = mark_cell_layers(sim, Vp, named_boundaries=['not left'])
    for cid in cells:
        cell = dolfin.Cell(mesh, cid)
        mp = cell.midpoint()[:]
        assert mp[0] > 0.9 or mp[1] < 0.1 or mp[1] > 0.9
    assert mpi_int_sum(len(cells)) == 20 * 1 + 18 * 2

    cells = mark_cell_layers(sim, Vp, named_boundaries=['left', 'not left'])
    assert mpi_int_sum(len(cells)) == 20 * 2 + 16 * 2


def mpi_int_sum(value):
    v = dolfin.MPI.sum(dolfin.MPI.comm_world, float(value))
    return int(v)
