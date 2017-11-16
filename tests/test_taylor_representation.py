import numpy
import dolfin
from ocellaris.utils import lagrange_to_taylor, taylor_to_lagrange
import pytest


@pytest.mark.parametrize("degree", [1, 2])
def test_taylor_projections_2D(degree):
    mesh = dolfin.UnitSquareMesh(4, 4)
    V = dolfin.FunctionSpace(mesh, 'DG', degree)
    
    # Setup Lagrange function with random data
    f1 = dolfin.Function(V)
    N = f1.vector().local_size()
    f1.vector().set_local(numpy.random.rand(N))
    
    # Convert to Taylor representation
    f2 = dolfin.Function(V)
    lagrange_to_taylor(f1,  f2)
    
    # Create a new Lagrange function from the Taylor representation
    f3 = dolfin.Function(V)
    taylor_to_lagrange(f2, f3)
    
    error = dolfin.errornorm(f1, f3, degree_rise=0)
    assert error < 1e-15


@pytest.mark.parametrize("degree", [1, 2])
def test_taylor_projections_3D(degree):
    if degree == 2:
        return pytest.skip('Skip while fixing degree=1')
    mesh = dolfin.UnitCubeMesh(2, 2, 2)
    V = dolfin.FunctionSpace(mesh, 'DG', degree)
    
    # Setup Lagrange function with random data
    f1 = dolfin.Function(V)
    N = f1.vector().local_size()
    f1.vector().set_local(numpy.random.rand(N))
    
    # Convert to Taylor representation
    f2 = dolfin.Function(V)
    lagrange_to_taylor(f1,  f2)
    
    # Create a new Lagrange function from the Taylor representation
    f3 = dolfin.Function(V)
    taylor_to_lagrange(f2, f3)
    
    from ocellaris.utils.taylor_basis import CACHE
    A1 = CACHE[('lagrange_to_taylor_matrices', degree)]
    A2 = CACHE[('taylor_to_lagrange_matrices', degree)]
    for A1i, A2i in zip(A1, A2):
        Ii = numpy.dot(A1i, A2i)
        print(Ii)
    print()
    print(A1[1], numpy.linalg.inv(A2[1]),
          A2[1], numpy.linalg.inv(A1[1]),
          sep='\n')
    
    error = dolfin.errornorm(f1, f3, degree_rise=0)
    assert error < 1e-15
