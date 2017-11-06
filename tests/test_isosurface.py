import dolfin
import numpy
from ocellaris import Simulation, setup_simulation
import pytest


ISO_INPUT = """
ocellaris:
    type: input
    version: 1.0
mesh:
    type: Rectangle
    Nx: 4
    Ny: 4
probes:
-   name: free_surface
    enabled: yes
    type: IsoSurface
    value: 0.5
    field: c
    custom_hook: MultiPhaseModelUpdated
multiphase_solver:
    type: BlendedAlgebraicVOF
    function_space_colour: DG
    polynomial_degree_colour: 0
solver: {type: AnalyticalSolution}
boundary_conditions: [{'name': 'all', 'selector': 'code', 'inside_code': 'on_boundary'}]
time: {dt: 1.0}
physical_properties: {nu0: 1.0, nu1: 1, rho0: 1, rho1: 1}
"""

@pytest.mark.parametrize("degree", [0, 1, 2])
def test_isoline_horizontal(degree):
    sim = Simulation()
    sim.input.read_yaml(yaml_string=ISO_INPUT)
    sim.input.set_value('multiphase_solver/polynomial_degree_colour', degree)
    setup_simulation(sim)
    probe = sim.probes['free_surface']

    # Initial value with sharp interface at x[1] == 0.5
    Vc = sim.data['Vc']
    c = sim.data['c']
    dm = Vc.dofmap()
    arr = c.vector().get_local()
    for cell in dolfin.cells(sim.data['mesh']):
        cell_value = 1 if cell.midpoint().y() < 0.5 else 0
        for dof in dm.cell_dofs(cell.index()):
            arr[dof] = cell_value
    c.vector().set_local(arr)
    c.vector().apply('insert')
    
    lines = probe.run(force_active=True)
    print('\nDegree:', degree, 'Vcdim:', Vc.dim())
    print(probe.name, probe.field_name, probe.value)
    print(len(lines))
    for x, y in lines:
        print('x', x, '\ny', y)
        assert all(abs(y - 0.5) < 1e-12)
        
        # Results should be in sorted order
        xdx = numpy.diff(x)
        assert all(xdx > 0) or all(xdx < 0)
    assert len(lines) == 1
