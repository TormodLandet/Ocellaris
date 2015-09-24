# encoding: utf8
from dolfin import TrialFunction, TestFunction, FacetNormal, Function
from dolfin import ds, dS, dot, as_vector
from dolfin import assemble_system, solve, cells, facets
import numpy


def convert_to_dgt(w, w_hat, interpolate=True):
    """
    Change the DGT function w_hat such that it takes the value w⋅n at exterior
    facets and the upwind value of w⋅n⁺ at interior facets if w_hat is scalar
    valued. If w_hat is vector valued then w_hat will equal the normal
    component of w at exterior facets and normal component of the upwind value
    of w at interior facets.  
    
    The conversion can be done via interpolation if w is a DG function. This
    is the fastest option. Otherwise projection can be used which involves a
    global matrix assembly and solve.
    
    The w function is assumed to be created as as_vector([u0, u1]) while the
    w_hat function is created from a VectorFunctionSpace.
    """
    Vdg = w[0].function_space()
    Vdgt = w_hat.function_space()
    
    assert Vdg.mesh().hash() == Vdgt.mesh().hash()
    assert Vdg.mesh().topology().dim() == 2
    assert Vdgt.ufl_element().family() == 'Discontinuous Lagrange Trace'
    assert w.ufl_shape == (2,)
    assert Vdg.ufl_element().degree() == Vdgt.ufl_element().degree()
    
    if interpolate:
        assert Vdg.ufl_element().family() == 'Discontinuous Lagrange'
        assert Vdg.ufl_element().degree() == Vdgt.ufl_element().degree() 
        _interpolate_to_dgt(w, w_hat)
    else:
        _project_to_dgt(w, w_hat)


def _project_to_dgt(w, w_hat):
    """
    This worker function performs projection into DGT from DG
    See the documentation of convert_to_dgt().
    """
    Vdgt = w_hat.function_space()
    mesh = Vdgt.mesh()
    
    u = TrialFunction(Vdgt)
    v = TestFunction(Vdgt)
    n = FacetNormal(mesh)
    
    # Should we produce scalar w⋅n or vector w⋅n n
    make_vectors =  Vdgt.num_sub_spaces() != 0
    
    if make_vectors:
        # Upwind flux
        wn = dot(w, n)
        wf_ds = wn*n
        w_normal_UW = (wn + abs(wn))/2*n
        wf_dS = w_normal_UW('+') + w_normal_UW('-')
    else:
        wn = dot(w, n)
        wf_ds = wn
        w_normal_UW = (wn + abs(wn))/2
        wf_dS = w_normal_UW('+') - w_normal_UW('-')
    
    a = dot(u, v)*ds
    L = dot(wf_ds, v)*ds
    for R in '+-':
        a += dot(u(R), v(R))*dS
        L += dot(wf_dS, v(R))*dS
    
    A, b = assemble_system(a, L)
    solve(A, w_hat.vector(), b)


def _interpolate_to_dgt(w, w_hat):
    """
    This worker function performs interpolation into DGT from DG
    See the documentation of convert_to_dgt().
    """
    Vdg = w[0].function_space()
    Vdgt = w_hat.function_space()
    mesh = Vdgt.mesh()
    D = mesh.topology().dim()
    k = Vdgt.ufl_element().degree()
    
    # Should we produce scalar w⋅n or vector w⋅n n
    make_vectors =  Vdgt.num_sub_spaces() != 0
    
    dofmap_dg = Vdg.dofmap()
    dofmap_dgt = Vdgt.dofmap()
    vec_dgs = [w[0].vector().array(), w[1].vector().array()]
    vec_dgt = w_hat.vector().array()
    vec_dgt[:] = 0
    
    # Mapping from DGT dofs to DG dofs on a cell
    if k == 0:
        mapping = [0, 0, 0]
    elif k == 1:
        mapping = [1, 2, 0, 2, 0, 1]
    elif k == 2:
        mapping = [1, 2, 3, 0, 2, 4, 0, 1, 5]
    
    # Perform local map for each facet in each cell
    tmp_w = numpy.zeros((k+1, D), float)
    for cell in cells(mesh):
        cell_dofs_dg = dofmap_dg.cell_dofs(cell.index())
        cell_dofs_dgt = dofmap_dgt.cell_dofs(cell.index())
        for i, facet in enumerate(facets(cell)):
            connected_cells = facet.entities(D)
            
            # Get values from DG vector
            for j in range(k+1):
                for d in range(D):
                    idx_dg = cell_dofs_dg[mapping[i*(k+1) + j]]
                    tmp_w[j,d] = vec_dgs[d][idx_dg]
            
            # The facet normal
            n = cell.normal(i)
            n = (n.x(), n.y())
            
            # Write into DGT vector
            for j in range(k+1):
                normal_vel = tmp_w[j,0]*n[0] + tmp_w[j,1]*n[1]
                
                # Check interior facets
                if len(connected_cells) > 1:
                    # Only upwind values will be written
                    if normal_vel < 1e-14:
                        continue
                    
                    # We store w⋅n⁺, so we swap the direction if we have n⁻
                    if not make_vectors and cell.index() != connected_cells[0]:
                        normal_vel *= -1
                
                if make_vectors:
                    for d in range(D):
                        idx_dgt = cell_dofs_dgt[i*(k+1) + j + d*(D + 1)*(k+1)]
                        vec_dgt[idx_dgt] += normal_vel*n[d]
                
                else:
                    idx_dgt = cell_dofs_dgt[i*(k+1) + j]
                    vec_dgt[idx_dgt] += normal_vel
                    
    w_hat.vector()[:] = vec_dgt


###############################################################################
## Tests

def test_interpolation_to_dgt_scalar(vector=False):
    from dolfin import UnitSquareMesh, VectorFunctionSpace, FunctionSpace
    import time
    
    N = 4
    k = 2
    mesh = UnitSquareMesh(N, N, 'right')
    
    if vector:
        V = VectorFunctionSpace(mesh, 'DGT', k)
    else:
        V = FunctionSpace(mesh, 'DGT', k)
    
    # Checkerboard initialization of known function w
    Vdg = FunctionSpace(mesh, 'DG', k)
    w = as_vector([Function(Vdg), Function(Vdg)]) 
    numpy.random.seed(42)
    for dof in range(Vdg.dim()):
        val = 1 + dof % 3 + numpy.random.rand()*0.1
        w[0].vector()[dof] = val
        w[1].vector()[dof] = val
    
    t1 = time.time()
    u_hat1 = Function(V)
    convert_to_dgt(w, u_hat1, interpolate=False)
    print 'Project takes %.3f seconds' % (time.time()-t1)
    
    t1 = time.time()
    u_hat2 = Function(V)
    convert_to_dgt(w, u_hat2, interpolate=True)
    print 'Interpolate takes %.3f seconds' % (time.time()-t1)
    
    rounded = lambda x: numpy.array([round(v, 8) for v in x.vector().array()][:10])
    print
    print rounded(u_hat1)
    print rounded(u_hat2)
    print
    
    norm = (u_hat1.vector().array()**2).sum()
    err = ((u_hat1.vector().array() - u_hat2.vector().array())**2).sum()
    print 'ERROR:', err, err/norm
    assert err/norm < 1e-20


def test_interpolation_to_dgt_vector():
    test_interpolation_to_dgt_scalar(True)


if __name__ == '__main__':
    for func in (test_interpolation_to_dgt_scalar,
                 test_interpolation_to_dgt_vector):
        print '#'*80
        print func.__name__
        print
        func()
        print
