from __future__ import division, print_function
import time, subprocess, os, sys
from math import log
import dolfin
from ocellaris import Simulation, setup_simulation, run_simulation


def run_and_calculate_error(N, dt, tmax, polydeg_u, polydeg_p, modifier=None):
    """
    Run Ocellaris and return L2 & H1 errors in the last time step
    """
    say(N, dt, tmax, polydeg_u, polydeg_p)

    # Setup and run simulation
    sim = Simulation()
    sim.input.read_yaml('taylor-green.inp')

    if sim.input.get_value('mesh/type') == 'Rectangle':
        # Use structured mesh
        sim.input.set_value('mesh/Nx', N)
        sim.input.set_value('mesh/Ny', N)
    else:
        # Create unstructured mesh with gmsh
        cmd1 = [
            'gmsh',
            '-string',
            'lc = %f;' % (2.0 / N),
            '-o',
            'taylor-green_%d.msh' % N,
            '-2',
            'taylor-green.geo',
        ]
        cmd2 = ['dolfin-convert', 'taylor-green_%d.msh' % N, 'taylor-green.xml']
        with open('/dev/null', 'w') as devnull:
            for cmd in (cmd1, cmd2):
                say(' '.join(cmd))
                subprocess.call(cmd, stdout=devnull, stderr=devnull)

    sim.input.set_value('time/dt', dt)
    sim.input.set_value('time/tmax', tmax)
    sim.input.set_value('solver/polynomial_degree_velocity', polydeg_u)
    sim.input.set_value('solver/polynomial_degree_pressure', polydeg_p)
    sim.input.set_value('output/stdout_enabled', False)

    if sim.input.get_value('solver/timestepping_method', 'BDF') == 'CN':
        sim.input.set_value(
            'initial_conditions/p/cpp_code',
            '-(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*(t+dt/2))/4',
        )

    # Turn off BDM
    # sim.input.set_value('solver/velocity_postprocessing', 'None')

    if modifier:
        modifier(sim)  # Running regression tests, modify some input params

    say(
        'Running with %s %s solver ...'
        % (
            sim.input.get_value('solver/type'),
            sim.input.get_value('solver/function_space_velocity', 'DG'),
        )
    )
    t1 = time.time()
    setup_simulation(sim)
    if 'Vcoupled' in sim.data:
        say('Num unknowns', sim.data['Vcoupled'].dim())
    run_simulation(sim)
    duration = time.time() - t1
    say('DONE')

    # Interpolate the analytical solution to the same function space
    Vu = sim.data['Vu']
    Vp = sim.data['Vp']
    vals = dict(
        t=sim.time,
        dt=sim.dt,
        nu=sim.input['physical_properties']['nu'],
        rho=sim.input['physical_properties']['rho'],
    )
    u0e = dolfin.Expression(
        sim.input.get_value('initial_conditions/up0/cpp_code'), degree=polydeg_u + 3, **vals
    )
    u1e = dolfin.Expression(
        sim.input.get_value('initial_conditions/up1/cpp_code'), degree=polydeg_u + 3, **vals
    )
    if sim.input.get_value('solver/timestepping_method', 'BDF') == 'CN':
        vals['t'] = sim.time - sim.dt
    pe = dolfin.Expression(
        sim.input.get_value('initial_conditions/p/cpp_code'), degree=polydeg_p + 3, **vals
    )

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

    say('Number of time steps:', sim.timestep)
    loglines = sim.log.get_full_log().split('\n')
    say('Num inner iterations:', sum(1 if 'iteration' in line else 0 for line in loglines))
    int_p = dolfin.assemble(sim.data['p'] * dolfin.dx)
    say('Number of mesh cells:', sim.data['mesh'].num_cells())
    say('p*dx', int_p)
    div_u_Vp = abs(dolfin.project(dolfin.div(sim.data['u']), Vp).vector().get_local()).max()
    say('div(u)|Vp', div_u_Vp)
    div_u_Vu = abs(dolfin.project(dolfin.div(sim.data['u']), Vu).vector().get_local()).max()
    say('div(u)|Vu', div_u_Vu)
    Vdg0 = dolfin.FunctionSpace(sim.data['mesh'], "DG", 0)
    div_u_DG0 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg0).vector().get_local()).max()
    say('div(u)|DG0', div_u_DG0)
    Vdg1 = dolfin.FunctionSpace(sim.data['mesh'], "DG", 1)
    div_u_DG1 = abs(dolfin.project(dolfin.div(sim.data['u']), Vdg1).vector().get_local()).max()
    say('div(u)|DG1', div_u_DG1)

    if 'u_mesh' in sim.data:
        Vmesh = sim.data['Vmesh']
        div_u_mesh_Vmesh = abs(
            dolfin.project(dolfin.div(sim.data['u_mesh']), Vmesh).vector().get_local()
        ).max()
        say('div(u_mesh)|V_mesh', div_u_mesh_Vmesh)
        div_u_mesh_DG0 = abs(
            dolfin.project(dolfin.div(sim.data['u_mesh']), Vdg0).vector().get_local()
        ).max()
        say('div(u_mesh)|DG0', div_u_mesh_DG0)
        div_u_mesh_DG1 = abs(
            dolfin.project(dolfin.div(sim.data['u_mesh']), Vdg1).vector().get_local()
        ).max()
        say('div(u_mesh)|DG1', div_u_mesh_DG1)

    if False:
        # Plot the results
        for fa, name in ((u0a, 'u0'), (u1a, 'u1'), (pa, 'p')):
            p1 = dolfin.plot(sim.data[name] - fa, title='%s_diff' % name, key='%s_diff' % name)
            p2 = dolfin.plot(fa, title=name + ' analytical', key=name)
            p1.write_png('%g_%g_%s_diff' % (N, dt, name))
            p2.write_png('%g_%g_%s' % (N, dt, name))
        dolfin.interactive()

    if N == 40 and False:
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
    for normname, selected in [('L2', slice(0, 3)), ('H1', slice(3, 6))]:
        say('======= ========== ========== ========== ===== ===== ===== ========= =====')
        say(' Discr.        Errors in %s norm         Convergence rates     Duration   ' % normname)
        say('------- -------------------------------- ----------------- ---------------')
        say('    %3s         u0         u1          p    u0    u1     p wallclock  rate' % restype)
        say('======= ========== ========== ========== ===== ===== ===== ========= =====')
        for i, idx in enumerate(indices):
            if idx not in results:
                break
            hmin, dt, duration = results[idx][-3:]
            eu0, eu1, ep = results[idx][selected]
            discr = hmin if restype == 'h' else dt
            say('%7.5f %10.2e %10.2e %10.2e' % (discr, eu0, eu1, ep), end=' ')
            if i > 0:
                prev_idx = indices[i - 1]
                prev_eu0, prev_eu1, prev_ep = results[prev_idx][selected]
                prev_hmin, prev_dt, prev_duration = results[prev_idx][-3:]
                prev_discr = prev_hmin if restype == 'h' else prev_dt
                fac = log(prev_discr / discr)
                say(
                    '%5.2f %5.2f %5.2f'
                    % (
                        log(prev_eu0 / eu0) / fac,
                        log(prev_eu1 / eu1) / fac,
                        log(prev_ep / ep) / fac,
                    ),
                    end=' ',
                )
                say(
                    '%9s %5.2f' % (seconds_as_string(duration), log(duration / prev_duration) / fac)
                )
            else:
                say('                  %9s' % seconds_as_string(duration))

        say('======= ========== ========== ========== ===== ===== ===== ========= =====')
        say()


def seconds_as_string(seconds):
    mins, secs = seconds // 60, seconds % 60
    if not mins:
        return '%4.1fs' % secs
    else:
        return '%2dm %4.1fs' % (mins, secs)


def say(*args, **kwargs):
    if dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        print(*args, **kwargs)


def run_convergence_space(N_list, modifier=None):
    dt = 0.01
    tmax = 1.0
    results = {}
    prev_N = None
    for N in N_list:
        say('Running N = %g with dt = %g' % (N, dt))
        results[N] = run_and_calculate_error(
            N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1, modifier=modifier
        )
        print_results(results, N_list, 'h')


def run_convergence_time(dt_list, modifier=None):
    N = 200
    tmax = 6.0
    results = {}
    for dt in dt_list:
        t1 = time.time()
        say('Running dt =', dt)
        results[dt] = run_and_calculate_error(
            N=N, dt=dt, tmax=tmax, polydeg_u=2, polydeg_p=1, modifier=modifier
        )
        print_results(results, dt_list, 'dt')


def make_quasi_3D(sim):
    """
    Quasi 3D - run the same 2D Taylor-Green, but on a 3D mesh with
    Neumann conditions in the z-direction so simthat u3 should allways
    be zero
    """
    say('Making Quasi 3D simulation')

    inp = sim.input

    # Change mesh
    Nx = inp.get_value('mesh/Nx')
    depth = 2 / Nx
    inp.set_value('mesh/Nz', 1)
    inp.set_value('mesh/endz', depth)
    inp.set_value('mesh/type', 'Box')

    # Fix gravity
    inp.set_value('physical_properties/g', [0, 0, 0])

    # Modify existing pure Dirichlet BC
    inp.get_value('boundary_conditions/0/u/cpp_code').append('0.0')  # Z-vel on boundary

    # Add new BC with pure Neumann BCs
    bc2 = {
        'name': 'back and front',
        'selector': 'code',
        'inside_code': 'on_boundary and (x[2] < 1e-5 or x[2] > %r - 1e-5)' % depth,
        'u0': {'type': 'ConstantGradient', 'value': 0},
        'u1': {'type': 'ConstantGradient', 'value': 0},
        'u2': {'type': 'ConstantValue', 'value': 0},
        'p': {'type': 'ConstantGradient', 'value': 0},
    }
    inp.get_value('boundary_conditions').append(bc2)

    inp.set_value('solver/coupled', {'solver': 'lu', 'lu_method': 'mumps'})


if __name__ == '__main__':
    mods = []

    def modifier(sim):
        for mod in mods:
            mod(sim)

    if '--Q3D' in sys.argv:
        mods.append(make_quasi_3D)

    if '--verbose' in sys.argv:
        mods.append(lambda sim: sim.input.set_value('output/stdout_enabled', True))

    run_convergence_space([8, 16, 24, 32, 40], modifier)
    # run_convergence_time([5e-1, 2.5e-1, 1.25e-1, 6.25e-2, 3.12e-2])
    # run_convergence_time([2, 1, 0.5, 0.25, 0.125])
    # dolfin.interactive()

