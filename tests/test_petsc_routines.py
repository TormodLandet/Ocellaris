import numpy
import dolfin
from dolfin_utils.test import skip_in_parallel, true_false_fixture
from ocellaris.utils import matmul, create_block_matrix


# Generate empty matrix using the block matrix routine 
use_block_matrix = true_false_fixture

# Use petc4py or PETSc C++ routines
use_cpp = true_false_fixture


def mk_mat(mesh_size=1, order=2, block=False, use_cpp=True):
    mesh = dolfin.UnitSquareMesh(mesh_size, mesh_size)
    V = dolfin.FunctionSpace(mesh, 'DG', order)
    
    if block:
        A = create_block_matrix(V, order*3, use_cpp=use_cpp)
    else:
        u, v = dolfin.TrialFunction(V), dolfin.TestFunction(V)    
        A = dolfin.assemble(u*v*dolfin.dx)
        A.zero()
    
    return A


@skip_in_parallel
def test_matmul(use_block_matrix, use_cpp):
    indices = numpy.array([1, 2, 4], numpy.intc)
    blockA = numpy.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]], float)
    blockB = numpy.identity(3, float)
    
    A = mk_mat(block=use_block_matrix, use_cpp=use_cpp)
    B = mk_mat(block=use_block_matrix, use_cpp=use_cpp)
    
    A.set(blockA, indices, indices)
    B.set(blockB, indices, indices)
    A.apply('insert')
    B.apply('insert')
    
    dolfin.parameters['linear_algebra_backend'] = 'PETSc'
    A = dolfin.as_backend_type(A)
    B = dolfin.as_backend_type(B)    
    print('A:\n', A.array())
    print('B:\n', B.array())
    
    C = matmul(A, B, use_cpp=use_cpp)
    print('C:\n', C.array())
    
    assert A.rank() == B.rank() == C.rank()
    assert A.size(0) == B.size(0) == C.size(0)
    assert A.size(1) == B.size(1) == C.size(1)
    
    Carr = C.array()
    Cnpy = numpy.dot(A.array(), B.array())
    assert (abs(Carr - Cnpy) < 1e-10).all()


if __name__ == '__main__':
    test_matmul(False)