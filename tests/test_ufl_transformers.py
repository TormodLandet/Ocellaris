import numpy
import ufl
from ufl.classes import Zero
import dolfin
from dolfin import UnitSquareMesh, FunctionSpace, VectorFunctionSpace, MixedElement
from dolfin import TestFunction, TestFunctions, TrialFunction, TrialFunctions
from dolfin import assemble, as_vector, Constant, dot, grad, dx
from ocellaris.utils import is_zero_ufl_expression, split_form_into_matrix
import pytest


def check_is_zero(expr, expected, verbose=True):
    val = is_zero_ufl_expression(expr, return_val=True)
    
    if verbose:
        print('is_zero_ufl_expression', val, expected, end=' ')
        print('   expr: %s    expr-repr: %r' % (expr, expr))
    
    val2 = 1 if val != 0 else 0
    assert val2 == expected


def test_is_zero_simple_scalar_expressions():
    mesh = UnitSquareMesh(4, 4)
    V = FunctionSpace(mesh, 'CG', 1)
    v = TestFunction(V)
    u = TrialFunction(V)
    
    check_is_zero(Zero()*u*v, 0)
    check_is_zero(Constant(0)*u*v, 1)
    check_is_zero(Zero()*v, 0)
    check_is_zero(Constant(0)*v, 1)
    check_is_zero(0*u*v, 0)
    check_is_zero(0*v, 0)
    check_is_zero(ufl.sin(0)*v, 0)
    check_is_zero(ufl.cos(0)*v, 1)


def test_is_zero_simple_vector_expressions():
    mesh = UnitSquareMesh(4, 4)
    V = FunctionSpace(mesh, 'CG', 1)
    v = TestFunction(V)
    u = TrialFunction(V)
    
    check_is_zero(dot(as_vector([Zero(), u]), as_vector([Zero(), v])), 1)
    check_is_zero(dot(as_vector([Zero(), u]), as_vector([v, Zero()])), 0)


def test_is_zero_with_nabla():
    mesh = UnitSquareMesh(4, 4)
    V1 = FunctionSpace(mesh, 'CG', 1)
    V2 = VectorFunctionSpace(mesh, 'CG', 1)
    v1 = TestFunction(V1)
    v2 = TrialFunction(V2)
    n = Constant([1, 1])
    nn = dolfin.outer(n, n)
    vv = dolfin.outer(v2, v2)
    
    check_is_zero(dot(n, grad(v1)), 1)
    check_is_zero(dolfin.div(v2), 1)
    check_is_zero(dot(dolfin.div(nn), n), 0)
    check_is_zero(dot(dolfin.div(vv), n), 1)
    check_is_zero(dolfin.inner(nn, grad(v2)), 1)


@pytest.mark.parametrize('shape', [(2, 2), (3, 3)])
def test_form_splitter_matrices(shape):
    mesh = UnitSquareMesh(dolfin.MPI.comm_self, 3, 3)
    Vu = FunctionSpace(mesh, 'DG', 2)
    Vp = FunctionSpace(mesh, 'DG', 1)

    def define_eq(u, v, p, q):
        "A simple Stokes-like coupled weak form"
        eq = Constant(2)*dot(u, v)*dx
        eq += dot(grad(p), v)*dx
        eq -= dot(Constant([1, 1]), v)*dx
        eq += dot(grad(q), u)*dx
        eq -= dot(Constant(0.3), q)*dx
        return eq
    
    if shape == (2, 2):
        eu = MixedElement([Vu.ufl_element(), Vu.ufl_element()])
        ew = MixedElement([eu, Vp.ufl_element()])
        W = FunctionSpace(mesh, ew)
        
        u, p = TrialFunctions(W)
        v, q = TestFunctions(W)
        u = as_vector([u[0], u[1]])
        v = as_vector([v[0], v[1]])
    
    elif shape == (3, 3):
        ew = MixedElement([Vu.ufl_element(), Vu.ufl_element(), Vp.ufl_element()])
        W = FunctionSpace(mesh, ew)
        
        u0, u1, p = TrialFunctions(W)
        v0, v1, q = TestFunctions(W)
        u = as_vector([u0, u1])
        v = as_vector([v0, v1])
        eq = define_eq(u, v, p, q)
    
    # Define coupled problem
    eq = define_eq(u, v, p, q)

    # Split the weak form into a saddle point block matrix system
    mat, vec = split_form_into_matrix(eq, W, W)

    # Check shape and that appropriate elements are none
    assert mat.shape == shape
    assert mat[-1,-1] is None
    for i in range(shape[0]-1):
        for j in range(shape[1]-1):
            if i != j:
                assert mat[i,j] is None
    
    # Check that the split matrices are identical with the
    # coupled matrix when reassembled into a single matrix
    compare_split_matrices(eq, mat, vec, W, W)


def compare_split_matrices(eq, mat, vec, Wv, Wu, eps=1e-14):
    # Assemble coupled system
    M, v = dolfin.system(eq)
    M = assemble(M).array()
    v = assemble(v).get_local()
    Ni, Nj = mat.shape
    
    def compute_diff(coupled, split):
        diff = abs(coupled - split).flatten().max()
        return diff
    
    # Rebuild coupled system from split parts
    M2 = numpy.zeros_like(M)
    v2 = numpy.zeros_like(v)
    for i in range(Ni):
        dm_Wi = Wv.sub(i).dofmap()
        
        if vec[i] is not None:
            data = assemble(vec[i]).get_local()
            dm_Vi = vec[i].arguments()[0].function_space().dofmap()
            for cell in dolfin.cells(Wv.mesh()):
                dofs_Wi = dm_Wi.cell_dofs(cell.index())
                dofs_Vi = dm_Vi.cell_dofs(cell.index())
                v2[dofs_Wi] = data[dofs_Vi]
        
        dofs_i = dm_Wi.dofs()
        diff = compute_diff(v[dofs_i], v2[dofs_i])
        print('Vector part %d has error %10.3e' % (i, diff),
              '<---' if diff > eps else '')

        for j in range(Nj):
            dm_Wj = Wv.sub(j).dofmap()
            if mat[i,j] is not None:
                data = assemble(mat[i,j]).array()
                dm_Vi = mat[i,j].arguments()[0].function_space().dofmap()
                dm_Vj = mat[i,j].arguments()[1].function_space().dofmap()
                for cell in dolfin.cells(Wv.mesh()):
                    dofs_Wi = dm_Wi.cell_dofs(cell.index())
                    dofs_Wj = dm_Wj.cell_dofs(cell.index())
                    dofs_Vi = dm_Vi.cell_dofs(cell.index())
                    dofs_Vj = dm_Vj.cell_dofs(cell.index())
                    W_idx = numpy.ix_(dofs_Wi, dofs_Wj)
                    V_idx = numpy.ix_(dofs_Vi, dofs_Vj)
                    M2[W_idx] = data[V_idx]
        
            dofs_j = dm_Wj.dofs()
            idx = numpy.ix_(dofs_i, dofs_j)
            diff = compute_diff(M[idx], M2[idx])
            print('Matrix part (%d, %d) has error %10.3e' % (i, j, diff),
                  '<---' if diff > eps else '')
    
    # Check that the original and rebuilt systems are identical
    assert compute_diff(M, M2) < eps
    assert compute_diff(v, v2) < eps
