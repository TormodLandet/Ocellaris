from dolfin import *
import numpy
import time


def project_to_dgt(w, w_hat):
    """
    Change the DGT function w_hat such that it takes the value w⋅n at
    exterior facets and the upwind value of w⋅n⁺ at interior facets
    
    The w function is vector valued and w_hat is scalar valued
    """
    Vdg = w.function_space()
    Vdgt = w_hat.function_space()
    mesh = Vdgt.mesh()
    D = mesh.topology().dim()
    N = Vdg.num_sub_spaces()
    k = Vdgt.ufl_element().degree()
    
    assert Vdg.ufl_element().family() == 'Discontinuous Lagrange'
    assert Vdgt.ufl_element().family() == 'Discontinuous Lagrange Trace'
    assert D == 2
    assert k in (0, 1, 2)
    assert N == 2
    assert Vdg.ufl_element().degree() == k
    
    un = TrialFunction(Vdgt)
    v = TestFunction(Vdgt)
    n = FacetNormal(mesh)
    
    # Upwind flux
    wn = dot(w, n)
    w_normal_UW = (wn + abs(wn))/2
    wf = w_normal_UW('+') - w_normal_UW('-')
    
    a = un*v*ds
    L = wn*v*ds
    for R in '+-':
        a += un(R)*v(R)*dS
        L += wf*v(R)*dS

    A, b = assemble_system(a, L)
    solve(A, w_hat.vector(), b)


def interpolate_to_dgt(w, w_hat):
    """
    Change the DGT function u_hat such that it takes the value w⋅n at
    exterior facets and the upwind value of w⋅n⁺ at interior facets
    
    The w function is vector valued and w_hat is scalar valued
    """
    Vdg = w.function_space()
    Vdgt = w_hat.function_space()
    mesh = Vdgt.mesh()
    D = mesh.topology().dim()
    N = Vdg.num_sub_spaces()
    k = Vdgt.ufl_element().degree()
    
    assert Vdg.ufl_element().family() == 'Discontinuous Lagrange'
    assert Vdgt.ufl_element().family() == 'Discontinuous Lagrange Trace'
    assert D == 2
    assert k in (0, 1, 2)
    assert N == 2
    assert Vdg.ufl_element().degree() == k
    
    dofmap_dg = Vdg.dofmap()
    dofmap_dgt = Vdgt.dofmap()
    vec_dg = w.vector().array()
    vec_dgt = w_hat.vector()
    vec_dgt[:] = 0
    
    # Mapping from DGT dofs to DG dofs on a cell
    if k == 0:
        mapping = [0, 0, 0]
    elif k == 1:
        mapping = [1, 2, 0, 2, 0, 1]
    elif k == 2:
        mapping = [1, 2, 3, 0, 2, 4, 0, 1, 5]
    
    tmp_u = numpy.zeros((k+1, D), float)
    for cell in cells(mesh):
        cell_dofs_dg = dofmap_dg.cell_dofs(cell.index())
        cell_dofs_dgt = dofmap_dgt.cell_dofs(cell.index())
        for i, facet in enumerate(facets(cell)):
            connected_cells = facet.entities(D)
            
            # Get values from DG vector
            for j in range(k+1):
                for d in range(N):
                    idx_dg = cell_dofs_dg[d*k*(D+1) + mapping[i*(k+1) + j]]
                    tmp_u[j,d] = vec_dg[idx_dg]
            
            # The facet normal
            n = cell.normal(i)
            
            # Write into DGT vector
            for j in range(k+1):
                normal_vel = tmp_u[j,0]*n.x() + tmp_u[j,1]*n.y()
                
                # Check interior facets
                if len(connected_cells) > 1:
                    # Only upwind values will be written
                    if normal_vel < 1e-14:
                        continue
                    # We store w⋅n⁺, so we swap the direction if we have n⁻
                    if cell.index() != connected_cells[0]:
                        normal_vel *= -1
                
                idx_dgt = cell_dofs_dgt[i*(k+1) + j]
                vec_dgt[idx_dgt] = normal_vel


def test_interpolation_to_dgt():
    N = 4
    k = 2
    mesh = UnitSquareMesh(N, N, 'right')
    V = FunctionSpace(mesh, 'DGT', k)
    
    # Checkerboard initialization of known function w
    Vdg = VectorFunctionSpace(mesh, 'DG', k)
    w = Function(Vdg)
    for cell in cells(mesh):
        cell_dofs_dg = Vdg.dofmap().cell_dofs(cell.index())
        for dof1, dof2 in zip(cell_dofs_dg[:3*k], cell_dofs_dg[3*k:]):
            vr = numpy.random.rand() * 0.1
            w.vector()[dof1] = 1 + cell.index()%2 + vr
            w.vector()[dof2] = 1 + cell.index()%2 + vr
    
    t1 = time.time()
    u_hat1 = Function(V)
    project_to_dgt(w, u_hat1)
    print 'Project takes %.3f seconds' % (time.time()-t1)
    
    t1 = time.time()
    u_hat2 = Function(V)
    interpolate_to_dgt(w, u_hat2)
    print 'Interpolate takes %.3f seconds' % (time.time()-t1)
    
    rounded = lambda x: numpy.array([round(v, 8) for v in x.vector().array()][:10])
    
    print rounded(u_hat1)
    print rounded(u_hat2)
    err = abs(u_hat1.vector().array() - u_hat2.vector().array()).sum()
    print err
    assert err < 1e-13
