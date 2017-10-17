import ufl
from ufl.classes import Zero
import dolfin
from dolfin import UnitSquareMesh, FunctionSpace, MixedElement
from dolfin import TestFunction, TestFunctions, TrialFunction, TrialFunctions
from dolfin import as_vector, Constant, dot, grad, dx
from ocellaris.utils import is_zero_ufl_expression, split_form_into_matrix


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


def test_form_splitter_coupled():
    mesh = UnitSquareMesh(1, 1)
    Vu = FunctionSpace(mesh, 'CG', 2)
    Vp = FunctionSpace(mesh, 'CG', 1)
    We = MixedElement([Vu.ufl_element(), Vu.ufl_element(), Vp.ufl_element()])
    W = FunctionSpace(mesh, We)
    
    u0, u1, p = TrialFunctions(W)
    v0, v1, q = TestFunctions(W)
    u = as_vector([u0, u1])
    v = as_vector([v0, v1])
    
    # A simple Stokes-like coupled weak form
    eq = Constant(2)*dot(u, v)*dx
    eq += dot(grad(p), v)*dx
    eq -= dot(Constant([1, 1]), v)*dx
    eq += dot(grad(q), u)*dx
    
    # Split the weak form into a 3x3 block matrix system
    mat, vec = split_form_into_matrix(eq, W, W)
    ((A00, A01, B0),
     (A10, A11, B1),
     (C0,   C1,  D)) = mat
    E0, E1, F = vec
    
    # Check that the blocks can be assembled
    for name in 'A00 A01 A10 A11 B0 B1 C0 C1 D E0 E1 F'.split():
        ufl_form = locals()[name]
        
        print('%s:' % name, end=' ')
        if ufl_form is not None:
            jited_form = dolfin.Form(ufl_form)
            linalg_obj = dolfin.assemble(jited_form)
            
            if linalg_obj.rank() == 2:
                print (linalg_obj.size(0), linalg_obj.size(1)),
            else:
                print (linalg_obj.size(),),
    print(ufl_form)
    
    # No coupling between velocity components
    assert A01 is None and A10 is None
    
    # Check that this is a saddle point system
    assert D is None
    
    # The RHS of the divergence of u is zero
    assert F is None
    
    # The rest of the blocks should be present
    assert None not in (A00, A11, B0, B1, C0, C1, E0, E1)


if __name__ == '__main__':
    import logging
    ffc_logger = logging.getLogger('FFC')
    ffc_logger.setLevel(logging.WARNING)
    
    test_is_zero_simple_scalar_expressions()
    test_is_zero_simple_vector_expressions()
    
    test_form_splitter_coupled()
