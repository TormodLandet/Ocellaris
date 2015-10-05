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
    sim.input.read_yaml('taylor-green.inp')
    
    #sim.input.set_value('solver/type', 'IPCS')
    
    if sim.input.get_value('mesh/type') == 'Rectangle':
        sim.input['mesh']['Nx'] = sim.input['mesh']['Ny'] = N
    else:
        # Create unstructured mesh with gmsh
        cmd1 = ['gmsh', '-string', 'lc = %f;' % (2.0/N),
                '-o', 'taylor-green_%d.msh' % N, '-2', 'taylor-green.geo']
        cmd2 = ['dolfin-convert', 'taylor-green_%d.msh' % N, 'taylor-green.xml']
        with open('/dev/null', 'w') as devnull:
            for cmd in (cmd1, cmd2):
                print ' '.join(cmd)
                subprocess.call(cmd, stdout=devnull, stderr=devnull)
        
    sim.input['time']['dt'] = dt
    sim.input['time']['tmax'] = tmax
    sim.input['solver']['polynomial_degree_velocity'] = polydeg_u
    sim.input['solver']['polynomial_degree_pressure'] = polydeg_p
    
    if sim.input.get_value('solver/timestepping_method') == 'CN':
        sim.input.set_value('initial_conditions/p/cpp_code', '-(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*(t+dt/2))/4') 
    
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
    
    # Interpolate the analytical solution to the same function space
    Vu = sim.data['Vu']
    Vp = sim.data['Vp']
    vals = dict(t=sim.time, dt=sim.dt, nu=sim.input['physical_properties']['nu0'])
    u0e = dolfin.Expression(sim.input.get_value('initial_conditions/up0/cpp_code'), **vals)
    u1e = dolfin.Expression(sim.input.get_value('initial_conditions/up1/cpp_code'), **vals)
    if sim.input.get_value('solver/timestepping_method') == 'CN':
        vals['t'] = sim.time - sim.dt
    pe  = dolfin.Expression(sim.input.get_value('initial_conditions/p/cpp_code'), **vals)
    
    #vals = dict(t=sim.time, nu=sim.input['physical_properties']['nu0'])
    #u0e = dolfin.Expression('-sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)', **vals)
    #u1e = dolfin.Expression(' sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)', **vals)
    #pe  = dolfin.Expression('-(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*t)/4', **vals)
    
    u0a = dolfin.project(u0e, Vu)
    u1a = dolfin.project(u1e, Vu)
    pa = dolfin.project(pe, Vp)
    
    # Calculate L2 errors
    err_u0 = calc_err(sim.data['u0'], u0a)
    err_u1 = calc_err(sim.data['u1'], u1a)
    err_p = calc_err(sim.data['p'], pa)
    
    # Calculate H1 errors
    err_u0_H1 = calc_err(sim.data['u0'], u0a, 'H1')
    err_u1_H1 = calc_err(sim.data['u1'], u1a, 'H1')
    err_p_H1 = calc_err(sim.data['p'], pa, 'H1')
    
    print 'Num iterations:', sim.timestep
    int_p = dolfin.assemble(sim.data['p']*dolfin.dx)
    print 'p*dx', int_p
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
    tmax = 1.0
    results = {}
    prev_N = None
    for N in N_list:
        print 'Running N = %g with dt = %g' % (N, dt)
        results[N] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1)
        print_results(results, N_list, 'h')


def run_convergence_time(dt_list):
    N = 200
    tmax = 6.0
    results = {}
    for dt in dt_list:
        t1 = time.time()
        print 'Running dt =', dt
        results[dt] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1)
        print_results(results, dt_list, 'dt')


run_convergence_space([8, 16, 24])
#run_convergence_time([5e-1, 2.5e-1, 1.25e-1, 6.25e-2, 3.12e-2])
#run_convergence_time([2, 1, 0.5, 0.25, 0.125])
