from __future__ import division, print_function
import time, subprocess, os
from math import log
import numpy
from matplotlib import pyplot
import dolfin
from ocellaris import Simulation, setup_simulation, run_simulation


def run_and_calculate_error(N, dt, tmax, polydeg_u, polydeg_p, last=False):
    """
    Run Ocellaris and return L2 & H1 errors in the last time step
    """
    say(N, dt, tmax, polydeg_u, polydeg_p)
    
    # Setup and run simulation
    timingtypes = [dolfin.TimingType_user, dolfin.TimingType_system, dolfin.TimingType_wall]
    dolfin.timings(dolfin.TimingClear_clear, timingtypes)
    sim = Simulation()
    sim.input.read_yaml('disc.inp')
    
    mesh_type = sim.input.get_value('mesh/type')
    if mesh_type == 'XML':
        # Create unstructured mesh with gmsh
        cmd1 = ['gmsh', '-string', 'lc = %f;' % (3.14/N),
                '-o', 'disc_%d.msh' % N, '-2', 'disc.geo']
        cmd2 = ['dolfin-convert', 'disc_%d.msh' % N, 'disc.xml']
        with open('/dev/null', 'w') as devnull:
            for cmd in (cmd1, cmd2):
                say(' '.join(cmd))
                subprocess.call(cmd, stdout=devnull, stderr=devnull)
    elif mesh_type == 'UnitDisc':
        sim.input.set_value('mesh/N', N//2)
    else:
        sim.input.set_value('mesh/Nx', N)
        sim.input.set_value('mesh/Ny', N)
    
    sim.input.set_value('time/dt', dt)
    sim.input.set_value('time/tmax', tmax)
    sim.input.set_value('solver/polynomial_degree_velocity', polydeg_u)
    sim.input.set_value('solver/polynomial_degree_pressure', polydeg_p)
    sim.input.set_value('output/stdout_enabled', False)
    
    say('Running with %s %s solver ...' % (sim.input.get_value('solver/type'), sim.input.get_value('solver/function_space_velocity')))
    t1 = time.time()
    setup_simulation(sim)
    run_simulation(sim)
    duration = time.time() - t1
    say('DONE')
    
    # Interpolate the analytical solution to the same function space
    Vu = sim.data['Vu']
    Vp = sim.data['Vp']
    Vr = sim.data['Vrho']
    polydeg_r = Vr.ufl_element().degree()
    vals = dict(t=sim.time, dt=sim.dt, Q=sim.input.get_value('user_code/constants/Q'))
    rho_e  = dolfin.Expression(sim.input.get_value('initial_conditions/rho_p/cpp_code'), degree=polydeg_r, **vals)
    u0e = dolfin.Expression(sim.input.get_value('initial_conditions/up0/cpp_code'), degree=polydeg_u, **vals)
    u1e = dolfin.Expression(sim.input.get_value('initial_conditions/up1/cpp_code'), degree=polydeg_u, **vals)
    pe  = dolfin.Expression(sim.input.get_value('initial_conditions/p/cpp_code'), degree=polydeg_p, **vals)
    
    rho_a = dolfin.project(rho_e, Vr)
    u0a = dolfin.project(u0e, Vu)
    u1a = dolfin.project(u1e, Vu)
    pa = dolfin.project(pe, Vp)
    
    mesh = sim.data['mesh']
    n = dolfin.FacetNormal(mesh)
    
    # Correct for possible non-zero average p
    int_p = dolfin.assemble(sim.data['p']*dolfin.dx)
    int_pa = dolfin.assemble(pa*dolfin.dx)
    vol = dolfin.assemble(dolfin.Constant(1.0)*dolfin.dx(domain=mesh))
    pa_avg = int_pa/vol
    sim.data['p'].vector()[:] += pa_avg
    
    # Calculate L2 errors
    err_rho = calc_err(sim.data['rho'], rho_a)
    err_u0 = calc_err(sim.data['u0'], u0a)
    err_u1 = calc_err(sim.data['u1'], u1a)
    err_p = calc_err(sim.data['p'], pa)
    
    # Calculate H1 errors
    err_rho_H1 = calc_err(sim.data['rho'], rho_a, 'H1')
    err_u0_H1 = calc_err(sim.data['u0'], u0a, 'H1')
    err_u1_H1 = calc_err(sim.data['u1'], u1a, 'H1')
    err_p_H1 = calc_err(sim.data['p'], pa, 'H1')
    
    reports = sim.reporting.timestep_xy_reports
    say('Num time steps:', sim.timestep)
    say('Num cells:', mesh.num_cells())
    say('Co_max:', numpy.max(reports['Co']))
    say('Pe_max:', numpy.max(reports['Pe']))
    say('rho_min went from %r to %r' % (reports['min(rho)'][0], reports['min(rho)'][-1]))
    say('rho_max went from %r to %r' % (reports['max(rho)'][0], reports['max(rho)'][-1]))
    m0, m1 = reports['mass'][0], reports['mass'][-1]
    say('mass error %.3e (%.3e)' % (m1 - m0, (m1 - m0)/m0))
    say('vel repr error %.3e' % dolfin.assemble(dolfin.dot(sim.data['u'], n)*dolfin.ds))
    say('p*dx', int_p)
    div_u_Vp = abs(dolfin.project(dolfin.div(sim.data['u']), Vp).vector().array()).max()
    say('div(u)|Vp', div_u_Vp)
    div_u_Vu = abs(dolfin.project(dolfin.div(sim.data['u']), Vu).vector().array()).max()
    say('div(u)|Vu', div_u_Vu)
    Vdg0 = dolfin.FunctionSpace(mesh, "DG", 0)
    div_u_DG0 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg0).vector().array()).max()
    say('div(u)|DG0', div_u_DG0)
    Vdg1 = dolfin.FunctionSpace(mesh, "DG", 1)
    div_u_DG1 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg1).vector().array()).max()
    say('div(u)|DG1', div_u_DG1)
    
    isoparam = mesh.ufl_coordinate_element().degree() > 1
    allways_plot = True
    if (last or allways_plot) and (not isoparam or sim.input.get_value('mesh/type') == 'UnitDisc'):
        # Plot the results
        for fa, name in ((u0a, 'u0'), (u1a, 'u1'), (pa, 'p'), (rho_a, 'rho')):
            fh = sim.data[name]
            if isoparam:
                # Bug in matplotlib plotting for isoparametric elements
                mesh2 = dolfin.UnitDiscMesh(dolfin.mpi_comm_world(), N//2, 1, 2)
                ue = fa.function_space().ufl_element()
                V2 = dolfin.FunctionSpace(mesh2, ue.family(), ue.degree())
                fa2, fh2 = dolfin.Function(V2), dolfin.Function(V2)
                fa2.vector().set_local(fa.vector().get_local())
                fh2.vector().set_local(fh.vector().get_local())
                fa, fh = fa2, fh2
            discr = '' # '%g_%g_' % (N, dt)
            plot(fa, name + ' analytical', '%s%s_1analytical' % (discr, name))
            plot(fh, name + ' numerical', '%s%s_2numerical' % (discr, name))
            plot(fh - fa, name + ' diff', '%s%s_3diff' % (discr, name))
    
    hmin = mesh.hmin()
    return err_rho, err_u0, err_u1, err_p, err_rho_H1, err_u0_H1, err_u1_H1, err_p_H1, hmin, dt, duration


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


def plot(func, title, filename):
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    p = dolfin.plot(func, title=title, backend='matplotlib')
    fig.colorbar(p)
    fig.savefig(filename+'.png')
    pyplot.close(fig)


def print_results(results, indices, restype):
    for normname, selected in [('L2', slice(0, 4)),
                               ('H1', slice(4, 8))]:
        say('======= ========== ========== ========== ========== ===== ===== ===== ===== ========= =====')
        say(' Discr.             Errors in %s norm                  Convergence rates        Duration   ' % normname)
        say('------- ------------------------------------------- ----------------------- ---------------')
        say('    %3s        rho         u0         u1          p   rho    u0    u1     p wallclock  rate' % restype)
        say('======= ========== ========== ========== ========== ===== ===== ===== ===== ========= =====')
        for i, idx in enumerate(indices):
            if idx not in results:
                break
            hmin, dt, duration = results[idx][-3:]
            erho, eu0, eu1, ep = results[idx][selected]
            discr = hmin if restype == 'h' else dt
            say('%7.5f %10.2e %10.2e %10.2e %10.2e' % (discr, erho, eu0, eu1, ep), end=' ')
            if i > 0:
                prev_idx = indices[i-1]
                prev_erho, prev_eu0, prev_eu1, prev_ep = results[prev_idx][selected]
                prev_hmin, prev_dt, prev_duration = results[prev_idx][-3:]
                prev_discr = prev_hmin if restype == 'h' else prev_dt
                fac = log(prev_discr/discr) 
                say('%5.2f %5.2f %5.2f %5.2f' % (log(prev_erho/erho)/fac, log(prev_eu0/eu0)/fac,
                                                 log(prev_eu1/eu1)/fac, log(prev_ep/ep)/fac), end=' ')
                say('%9s %5.2f' % (seconds_as_string(duration), log(duration/prev_duration)/fac))
            else:
                say('                        %9s' % seconds_as_string(duration))
        
        say('======= ========== ========== ========== ========== ===== ===== ===== ===== ========= =====')
        say()


def seconds_as_string(seconds):
    mins, secs = seconds//60, seconds%60
    if not mins:
        return '%4.1fs' % secs
    else:
        return '%2dm %4.1fs' % (mins, secs)
    

def say(*args, **kwargs):
    if dolfin.MPI.rank(dolfin.mpi_comm_world()) == 0:
        print(*args, **kwargs)


def run_convergence_space(N_list):
    dt = 0.002
    tmax = numpy.pi/4
    results = {}
    prev_N = None
    for N in N_list:
        say('Running N = %g with dt = %g' % (N, dt))
        results[N] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1,
                                             last=(N == N_list[-1]))
        print_results(results, N_list, 'h')


def run_convergence_time(dt_list):
    N = 200
    tmax = 10.0
    results = {}
    for dt in dt_list:
        t1 = time.time()
        say('Running dt =', dt)
        results[dt] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1,
                                              last=(dt == dt_list[-1]))
        print_results(results, dt_list, 'dt')


run_convergence_space([4, 8, 16])#, 24])
#run_convergence_time([5e-1, 2.5e-1, 1.25e-1, 6.25e-2, 3.12e-2])
#run_convergence_time([2, 1, 0.5, 0.25, 0.125])
#dolfin.interactive()