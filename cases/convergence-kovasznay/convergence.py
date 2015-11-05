from __future__ import division
import time, subprocess, os
from math import log
import dolfin
from ocellaris import Simulation, run_simulation


def run_and_calculate_error(N, dt, tmax, polydeg_u, polydeg_p):
    """
    Run Ocellaris and return L2 & H1 errors in the last time step
    """
    print N, dt, tmax, polydeg_u, polydeg_p
    
    # Setup and run simulation
    sim = Simulation()
    sim.input.read_yaml('kovasznay.inp')
    
    sim.input.set_value('mesh/Nx', N)
    sim.input.set_value('mesh/Ny', N)
    sim.input.set_value('time/dt', dt)
    sim.input.set_value('time/tmax', tmax)
    sim.input.set_value('solver/polynomial_degree_velocity', polydeg_u)
    sim.input.set_value('solver/polynomial_degree_pressure', polydeg_p)
    sim.input.set_value('output/ocellaris_log_level', 'warning')
    
    print 'Running ...'
    try:
        t1 = time.time()
        run_simulation(sim)
        duration = time.time() - t1
    except KeyboardInterrupt:
        raise
    except BaseException as e:
        import traceback
        traceback.print_exc()
        return [1e10]*6 + [1, dt, time.time()-t1]
    print 'DONE'
    tmax_warning = ' <------ NON CONVERGENCE!!' if sim.time >= tmax else ''
    
    # Interpolate the analytical solution to the same function space
    Vu = sim.data['Vu']
    Vp = sim.data['Vp']
    u0e = dolfin.Expression(sim.input.get_value('boundary_conditions/0/u/cpp_code/0'), degree=polydeg_u)
    u1e = dolfin.Expression(sim.input.get_value('boundary_conditions/0/u/cpp_code/1'), degree=polydeg_u)
    pe  = dolfin.Expression('-0.5*exp(-0.9637405441957689*2*x[0])', degree=polydeg_p)
    u0a = dolfin.project(u0e, Vu)
    u1a = dolfin.project(u1e, Vu)
    pa = dolfin.project(pe, Vp)
    
    # Correct pa
    int_pa = dolfin.assemble(pa*dolfin.dx)
    vol = dolfin.assemble(dolfin.Constant(1.0)*dolfin.dx(domain=Vp.mesh()))
    pa.vector()[:] -= int_pa/vol
    
    # Calculate L2 errors
    err_u0 = calc_err(sim.data['u0'], u0a)
    err_u1 = calc_err(sim.data['u1'], u1a)
    err_p = calc_err(sim.data['p'], pa)
    
    # Calculate H1 errors
    err_u0_H1 = calc_err(sim.data['u0'], u0a, 'H1')
    err_u1_H1 = calc_err(sim.data['u1'], u1a, 'H1')
    err_p_H1 = calc_err(sim.data['p'], pa, 'H1')
    
    print 'Number of time steps:', sim.timestep, tmax_warning
    print 'max(ui_new-ui_prev)', sim.reporting.get_report('max(ui_new-ui_prev)')[1][-1]
    int_p = dolfin.assemble(sim.data['p']*dolfin.dx)
    print 'p*dx', int_p
    print 'pa*dx', dolfin.assemble(pa*dolfin.dx(domain=Vp.mesh()))
    div_u_Vp = abs(dolfin.project(dolfin.div(sim.data['u']), Vp).vector().array()).max()
    print 'div(u)|Vp', div_u_Vp
    div_u_Vu = abs(dolfin.project(dolfin.div(sim.data['u']), Vu).vector().array()).max()
    print 'div(u)|Vu', div_u_Vu
    Vdg0 = dolfin.FunctionSpace(sim.data['mesh'], "DG", 0)
    div_u_DG0 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg0).vector().array()).max()
    print 'div(u)|DG0', div_u_DG0
    Vdg1 = dolfin.FunctionSpace(sim.data['mesh'], "DG", 1)
    div_u_DG1 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg1).vector().array()).max()
    print 'div(u)|DG1', div_u_DG1
    
    if False:
        # Plot the results
        for fa, name in ((u0a, 'u0'), (u1a, 'u1'), (pa, 'p')): 
            p1 = dolfin.plot(sim.data[name] - fa, title='%s_diff' % name, key='%s_diff' % name)
            p2 = dolfin.plot(fa, title=name+' analytical', key=name)
            p1.write_png('%g_%g_%s_diff' % (N, dt, name))
            p2.write_png('%g_%g_%s' % (N, dt, name))
        dolfin.interactive()
        
    if N == 24 and False:
        dolfin.plot(sim.data['u0'], title='u0')
        dolfin.plot(sim.data['u1'], title='u1')
        dolfin.plot(sim.data['p'], title='p')
        dolfin.plot(u0a, title='u0a')
        dolfin.plot(u1a, title='u1a')
        dolfin.plot(pa, title='pa')
        plot_err(sim.data['u0'], u0a, 'u0a - u0')
        plot_err(sim.data['u1'], u1a, 'u1a - u1')
        plot_err(sim.data['p'], pa, 'pa - p')
    
    hmin = sim.data['mesh'].hmin()
    return err_u0, err_u1, err_p, err_u0_H1, err_u1_H1, err_p_H1, hmin, dt, duration


def calc_err(f_num, f_ana, normtype='l2'):
    """
    Calculate scaled L2 error
    """
    f_err = dolfin.Function(f_num.function_space())
    f_err.vector()[:] = f_ana.vector()[:] - f_num.vector()[:]
    if normtype == 'l2':
        return dolfin.norm(f_err) / dolfin.norm(f_ana)
    else:
        return dolfin.norm(f_err, normtype)


def plot_err(f_num, f_ana, title):
    f_err = dolfin.Function(f_num.function_space())
    f_err.vector()[:] = f_ana.vector()[:] - f_num.vector()[:]
    dolfin.plot(f_err, title=title)


def print_results(results, indices, restype):
    for normname, selected in [('L2', slice(0, 3)),
                               ('H1', slice(3, 6))]:
        print '======= ========== ========== ========== ===== ===== ===== ========= ====='
        print ' Discr.        Errors in %s norm         Convergence rates     Duration   ' % normname
        print '------- -------------------------------- ----------------- ---------------'
        print '    %3s         u0         u1          p    u0    u1     p wallclock  rate' % restype
        print '======= ========== ========== ========== ===== ===== ===== ========= ====='
        for i, idx in enumerate(indices):
            if idx not in results:
                break
            hmin, dt, duration = results[idx][-3:]
            eu0, eu1, ep = results[idx][selected]
            discr = hmin if restype == 'h' else dt
            print '%7.5f %10.2e %10.2e %10.2e' % (discr, eu0, eu1, ep),
            if i > 0:
                prev_idx = indices[i-1]
                prev_eu0, prev_eu1, prev_ep = results[prev_idx][selected]
                prev_hmin, prev_dt, prev_duration = results[prev_idx][-3:]
                prev_discr = prev_hmin if restype == 'h' else prev_dt
                fac = log(prev_discr/discr) 
                print '%5.2f %5.2f %5.2f' % (log(prev_eu0/eu0)/fac, log(prev_eu1/eu1)/fac, log(prev_ep/ep)/fac),
                print '%9s %5.2f' % (seconds_as_string(duration), log(duration/prev_duration)/fac)
            else:
                print '                  %9s' % seconds_as_string(duration)
        
        print '======= ========== ========== ========== ===== ===== ===== ========= ====='
        print


def seconds_as_string(seconds):
    mins, secs = seconds//60, seconds%60
    if not mins:
        return '%4.1fs' % secs
    else:
        return '%2dm %4.1fs' % (mins, secs)


def run_convergence_space(N_list):
    dt = 0.01
    tmax = 3.0
    results = {}
    prev_N = None
    for N in N_list:
        print 'Running N = %g with dt = %g' % (N, dt)
        results[N] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1)
        print_results(results, N_list, 'h')

run_convergence_space([8, 16, 24])
#dolfin.interactive()