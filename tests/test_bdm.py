import dolfin
from ocellaris import Simulation
from ocellaris.solver_parts import VelocityBDMProjection
import pytest


@pytest.mark.parametrize("gdim", [2, 3])
def test_bdm_identity_transform(gdim):
    N = 5
    k = 2

    # Create mesh and define a divergence free field
    if gdim == 2:
        mesh = dolfin.UnitSquareMesh(N, N)
        field = ['-sin(pi*x[1])*cos(pi*x[0])', 'sin(pi*x[0])*cos(pi*x[1])']
    elif gdim == 3:
        mesh = dolfin.UnitCubeMesh(N, N, N)
        field = [
            '-sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[2])',
            'sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[2])',
            '-cos(pi*x[2])*cos(pi*(x[0]-x[1]))',
        ]

    V = dolfin.FunctionSpace(mesh, 'DG', k)
    u = dolfin.as_vector([dolfin.Function(V) for _ in range(gdim)])
    a = dolfin.as_vector([dolfin.Function(V) for _ in range(gdim)])

    # Make a divergence free field
    for i, cpp in enumerate(field):
        ei = dolfin.Expression(cpp, degree=k)
        u[i].interpolate(ei)
        a[i].assign(u[i])

    # Run the projection into BDM and then assign this to u
    bdm = VelocityBDMProjection(simulation=Simulation(), w=u, use_bcs=False)
    bdm.run()

    # Check the changes made
    for ui, ai in zip(u, a):
        error = dolfin.errornorm(ai, ui, degree_rise=0)
        assert error < 1e-15
