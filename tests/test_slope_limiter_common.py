import dolfin
import pytest


@pytest.mark.parametrize('D', (2, 3))
def test_cell_midpoints(D):
    from ocellaris.solver_parts.slope_limiter.limiter_cpp_utils import SlopeLimiterInput 
    if D == 2:
        mesh = dolfin.UnitSquareMesh(4, 4)
    else:
        mesh = dolfin.UnitCubeMesh(2, 2, 2)
    
    Vx = dolfin.FunctionSpace(mesh, 'DG', 2)
    V0 = dolfin.FunctionSpace(mesh, 'DG', 0)
    
    py_inp = SlopeLimiterInput(mesh, Vx, V0)
    cpp_inp = py_inp.cpp_obj 
    
    all_ok = True
    for cell in dolfin.cells(mesh):
        cid = cell.index()
        mp = cell.midpoint()
        cpp_mp = cpp_inp.cell_midpoints[cid]
        for d in range(D):
            ok = dolfin.near(mp[d], cpp_mp[d])
            if not ok:
                print('%3d %d - %10.3e %10.3e' % (cid, d, mp[d], cpp_mp[d]), '<--- ERROR' if not ok else '')
                all_ok = False
    print(cid)
    
    assert all_ok
