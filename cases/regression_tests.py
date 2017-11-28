"""
This defines a set of small and relatively fast regression tests
which can be run by use of pytest. Run

  python3 -m pytest -vs cases/regresion_tests.py

to execute all regression tests
"""
import os, sys
import pytest


CASEDIR = os.path.abspath(os.path.dirname(__file__))
RESNAMES = 'err_u0 err_u1 err_p err_u0_H1 err_u1_H1 err_p_H1 hmin dt duration'.split()
SPEEDFOCUS = True


@pytest.fixture(params=['Coupled', 'SIMPLE'])
def solver_type(request):
    return request.param


def test_taylor_green(solver_type, monkeypatch):
    runner = regression_setup(monkeypatch, 'convergence-taylor-green')
    N = 8
    dt = 0.05
    tmax = 1
    polydeg_u = 2
    polydeg_p = 1
    limits = {'Coupled': [0.015, 0.015, 0.05, 0.28, 0.28, 0.18],
              'SIMPLE':  [0.017, 0.017, 0.05, 0.28, 0.28, 0.18]}
    
    def modifier(sim):
        sim.input.set_value('solver/type', solver_type)
        shared_regression_modifier(sim)
    
    # Run the Taylor-Green test case on a course mesh
    res = runner(N, dt, tmax, polydeg_u, polydeg_p, modifier)
    assert check_results(res, limits[solver_type])


def test_kovasznay(solver_type, monkeypatch):
    runner = regression_setup(monkeypatch, 'convergence-kovasznay')
    N = 8
    dt = 0.01
    tmax = 3.0
    polydeg_u = 2
    polydeg_p = 1
    limits = {'Coupled': [0.00037, 0.007, 0.02, 0.016, 0.023, 0.07],
              'SIMPLE':  [0.00500, 0.084, 0.08, 0.075, 0.130, 0.12]}
    
    def modifier(sim):
        sim.input.set_value('solver/type', solver_type)
        sim.input.set_value('solver/steady_velocity_stopping_criterion', 1e-5)
        shared_regression_modifier(sim)
    
    # Run the Taylor-Green test case on a course mesh 
    res = runner(N, dt, tmax, polydeg_u, polydeg_p, modifier)
    assert check_results(res, limits[solver_type])


####################################################################################
# Helpers

def regression_setup(monkeypatch, casename):
    casedir = os.path.join(CASEDIR, casename)
    assert os.path.isdir(casedir)
    monkeypatch.chdir(casedir)
    monkeypatch.syspath_prepend(casedir)
    from convergence import run_and_calculate_error
    sys.modules.pop('convergence')
    return run_and_calculate_error


def shared_regression_modifier(sim):
    sim.input.set_value('output/save_restart_file_at_end', False)
    sim.input.set_value('output/stdout_enabled', False)
    sim.input.set_value('output/xdmf_write_interval', 0)


def check_results(res, limits):
    OK = True
    for i, result in enumerate(res):
        name = RESNAMES[i]
        if i >= len(limits):
            print('%-10s: %g' % (name, result))
            continue
        
        limit = limits[i]
        if result > limit:
            OK = False
            print('%-10s: %g > %g  <--- ERROR' % (name, result, limit))
        else:
            print('%-10s: %g < %g' % (name, result, limit))
    
    return OK


