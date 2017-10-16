import os
import pytest


CASEDIR = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture(params=['Coupled', 'SIMPLE'])
def solver_type(request):
    return request.param


def test_taylor_green(solver_type, monkeypatch):
    monkeypatch.chdir(os.path.join(CASEDIR, 'convergence-taylor-green'))
    from convergence import run_and_calculate_error
    
    N = 8
    dt = 0.05
    tmax = 1
    polydeg_u = 2
    polydeg_p = 1
    
    # Run the Taylor-Green test case on a course mesh 
    err_u0, err_u1, err_p, err_u0_H1, err_u1_H1, err_p_H1, hmin, dt, duration \
      = run_and_calculate_error(N, dt, tmax, polydeg_u, polydeg_p, solver_type=solver_type)
    print(err_u0, err_u1, err_p, err_u0_H1, err_u1_H1, err_p_H1, hmin, dt, duration)
    
    assert err_u0 < 0.017
    assert err_u1 < 0.017
    assert err_p < 0.23
    assert err_u0_H1 < 0.31
    assert err_u1_H1 < 0.31
    assert err_p_H1 < 0.47
