from __future__ import division, print_function
import time, subprocess, os
from math import log
import numpy
from matplotlib import pyplot
import dolfin
from ocellaris import Simulation, setup_simulation, run_simulation


ISROOT = dolfin.MPI.rank(dolfin.MPI.comm_world) == 0


def run_and_calculate_error(N, dt, tmax, polydeg_rho, last=False):
    """
    Run Ocellaris and return L2 & H1 errors in the last time step
    """
    say(N, dt, tmax, polydeg_rho)
    
    # Setup and run simulation
    sim = Simulation()
    sim.input.read_yaml('transport.inp')
    
    mesh_type = sim.input.get_value('mesh/type')
    if mesh_type == 'XML':
        # Create unstructured mesh with gmsh
        cmd1 = ['gmsh', '-string', 'lc = %f;' % (3.14/N),
                '-o', 'disc_%d.msh' % N, '-2',
                '../convergence-variable-density-disk/disc.geo']
        cmd2 = ['dolfin-convert', 'disc_%d.msh' % N, 'disc.xml']
        with open('/dev/null', 'w') as devnull:
            for cmd in (cmd1, cmd2):
                say(' '.join(cmd))
                if ISROOT:
                    subprocess.call(cmd, stdout=devnull, stderr=devnull)
    elif mesh_type == 'UnitDisc':
        sim.input.set_value('mesh/N', N//2)
    else:
        sim.input.set_value('mesh/Nx', N)
        sim.input.set_value('mesh/Ny', N)
    
    sim.input.set_value('time/dt', dt)
    sim.input.set_value('time/tmax', tmax)
    sim.input.set_value('multiphase_solver/polynomial_degree_rho', polydeg_rho)
    sim.input.set_value('output/stdout_enabled', False)
    
    say('Running with multiphase solver %s ...' % (sim.input.get_value('multiphase_solver/type')))
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
    vals = dict(t=sim.time, dt=sim.dt)
    rho_e  = dolfin.Expression(sim.input.get_value('initial_conditions/rho_p/cpp_code'), degree=polydeg_r, **vals)
    rho_a = dolfin.project(rho_e, Vr)
    
    rho_e.t = 0
    rho_0 = dolfin.project(rho_e, Vr)
    
    # Calculate L2 errors
    err_rho = calc_err(sim.data['rho'], rho_a)
    
    # Calculate H1 errors
    err_rho_H1 = calc_err(sim.data['rho'], rho_a, 'H1')
    
    mesh = sim.data['mesh']
    n = dolfin.FacetNormal(mesh)
    
    reports = sim.reporting.timestep_xy_reports
    say('Num time steps:', sim.timestep)
    say('Num cells:', mesh.num_cells())
    say('Co_max:', numpy.max(reports['Co']))
    say('rho_min went from %r to %r' % (reports['min(rho)'][0], reports['min(rho)'][-1]))
    say('rho_max went from %r to %r' % (reports['max(rho)'][0], reports['max(rho)'][-1]))
    m0, m1 = reports['mass'][0], reports['mass'][-1]
    say('mass error %.3e (%.3e)' % (m1 - m0, (m1 - m0)/m0))
    say('vel compat error %.3e' % dolfin.assemble(dolfin.dot(sim.data['u'], n)*dolfin.ds))
    int_p = dolfin.assemble(sim.data['p']*dolfin.dx)
    say('p*dx', int_p)
    div_u_Vp = abs(dolfin.project(dolfin.div(sim.data['u']), Vp).vector().get_local()).max()
    say('div(u)|Vp', div_u_Vp)
    div_u_Vu = abs(dolfin.project(dolfin.div(sim.data['u']), Vu).vector().get_local()).max()
    say('div(u)|Vu', div_u_Vu)
    Vdg0 = dolfin.FunctionSpace(mesh, "DG", 0)
    div_u_DG0 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg0).vector().get_local()).max()
    say('div(u)|DG0', div_u_DG0)
    Vdg1 = dolfin.FunctionSpace(mesh, "DG", 1)
    div_u_DG1 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg1).vector().get_local()).max()
    say('div(u)|DG1', div_u_DG1)
    
    isoparam = mesh.ufl_coordinate_element().degree() > 1
    if last and (not isoparam or sim.input.get_value('mesh/type') == 'UnitDisc'):
        # Plot the results
        for fa, name in ((rho_a, 'rho'),):
            fh = sim.data[name]
            if isoparam:
                # Bug in matplotlib plotting for isoparametric elements
                mesh2 = dolfin.UnitDiscMesh(dolfin.MPI.comm_world, N//2, 1, 2)
                ue = fa.function_space().ufl_element()
                V2 = dolfin.FunctionSpace(mesh2, ue.family(), ue.degree())
                fa2, fh2 = dolfin.Function(V2), dolfin.Function(V2)
                fa2.vector().set_local(fa.vector().get_local())
                fh2.vector().set_local(fh.vector().get_local())
                fa, fh = fa2, fh2
            plot(fh - fa, name + ' diff', '%g_%g_%s_diff' % (N, dt, name))
            plot(fa, name + ' analytical', '%g_%g_%s_analytical' % (N, dt, name))
            plot(fh, name + ' numerical', '%g_%g_%s_numerical' % (N, dt, name))
            plot(rho_0, name + ' initial', '%g_%g_%s_initial' % (N, dt, name))
    
    hmin = mesh.hmin()
    return err_rho, err_rho_H1, hmin, dt, duration


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
    for normname, selected in [('L2', slice(0, 1)),
                               ('H1', slice(1, 2))]:
        say('======= =================== ================== ========= =====')
        say(' Discr.  Errors in %s norm   Convergence rate        Duration   ' % normname)
        say('------- ------------------- ------------------ ---------------')
        say('    %3s        rho                  rho        wallclock  rate' % restype)
        say('======= =================== ================== ========= =====')
        for i, idx in enumerate(indices):
            if idx not in results:
                break
            hmin, dt, duration = results[idx][-3:]
            erho, = results[idx][selected]
            discr = hmin if restype == 'h' else dt
            say('%7.5f      %10.2e    ' % (discr, erho), end=' ')
            if i > 0:
                prev_idx = indices[i-1]
                prev_erho, = results[prev_idx][selected]
                prev_hmin, prev_dt, prev_duration = results[prev_idx][-3:]
                prev_discr = prev_hmin if restype == 'h' else prev_dt
                fac = log(prev_discr/discr) 
                say('      %5.2f       ' % (log(prev_erho/erho)/fac), end=' ')
                say('%9s %5.2f' % (seconds_as_string(duration), log(duration/prev_duration)/fac))
            else:
                say('                   %9s' % seconds_as_string(duration))
        
        say('======= =================== ================== ========= =====')
        say()


def seconds_as_string(seconds):
    mins, secs = seconds//60, seconds%60
    if not mins:
        return '%4.1fs' % secs
    else:
        return '%2dm %4.1fs' % (mins, secs)
    

def say(*args, **kwargs):
    if ISROOT:
        print(*args, **kwargs)


def run_convergence_space(N_list):
    dt = 0.01
    P = 1
    tmax = numpy.pi/4
    results = {}
    prev_N = None
    for N in N_list:
        dx = 2/N/5
        Co = 1/(2*P+1)
        dt = Co*dx/2**0.5
        say('Running N = %g with dt = %g' % (N, dt))
        results[N] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_rho=P,
                                             last=(N == N_list[-1]))
        print_results(results, N_list, 'h')


def run_convergence_time(dt_list):
    N = 100
    tmax = 10.0
    results = {}
    for dt in dt_list:
        t1 = time.time()
        say('Running dt =', dt)
        results[dt] = run_and_calculate_error(N=N, dt=dt, tmax=tmax, polydeg_rho=3,
                                              last=(dt == dt_list[-1]))
        print_results(results, dt_list, 'dt')


run_convergence_space([4, 8, 16, 24])#, 48, 96])
#run_convergence_time([5e-1, 2.5e-1, 1.25e-1])#, 6.25e-2, 3.12e-2])
#run_convergence_time([2, 1, 0.5, 0.25, 0.125])
#dolfin.interactive()