from __future__ import division, print_function
import sys
import numpy
from matplotlib import pyplot
from matplotlib.tri import Triangulation
import dolfin as df
from ocellaris import Simulation
from plot_streamlines import load_simulation, StreamFunction


H1 = 1
H2 = 0.9423
L1, L2 = 5, 20
N_levels = 30


def plot_backward_facing_step(res):
    u = df.as_vector([res['u0'], res['u1']])
    V = u[0].function_space()
    mesh = V.mesh()

    ############################################################################
    # Get some meta data
    sim = Simulation()
    sim.input.read_yaml(yaml_string=res['input'])
    Re = sim.input.get_value('user_code/constants/Re', required_type='float')
    time = res['time']

    ############################################################################
    # Find the length X1 of primary recirculating bubble

    x1_ypos = -0.99 * H2
    xf, _, _, f = get_probe_values(u[0], [0, x1_ypos, 0], [L2, x1_ypos, 0], 200)
    x1_c0 = get_zero_crossings(xf, f)
    print(x1_c0)

    ############################################################################
    # Find the inflow profile

    x_inflow = -1.5 * H1
    _, y_inflow, _, u0_inflow = get_probe_values(
        u[0], [x_inflow, 0, 0], [x_inflow, H1, 0], 20
    )

    ############################################################################
    # Refine the mesh and get a triangluated stream function
    for _ in range(2):
        old = mesh, V
        mesh = df.refine(mesh)
        V = df.FunctionSpace(mesh, V.ufl_element())
        df.parameters['allow_extrapolation'] = True
        try:
            u0 = df.interpolate(u[0], V)
            u1 = df.interpolate(u[1], V)
            u = df.as_vector([u0, u1])
        except Exception:
            mesh, V = old
            u0, u1 = u[0], u[1]
            print('Refining failed ------------------------------')
            break
        df.parameters['allow_extrapolation'] = False

    sf = StreamFunction(u, degree=V.ufl_element().degree())
    sf.compute()

    levels = get_probe_values(
        sf.psi, [1.5 * H2, -0.5 * H2, 0], [H2, H1 * 0.99, 0], N_levels
    )[-1]
    levels2 = []
    for xi, _idx, _upcrossing in x1_c0:
        # Add stream function values at upcrossing locations
        if xi - 0.3 > 0:
            levels2.append(get_probe_value(sf.psi, (xi - 0.3, x1_ypos, 0)))
        if xi + 0.1 < L2:
            levels2.append(get_probe_value(sf.psi, (xi + 0.1, x1_ypos, 0)))
    if levels2:
        levels = numpy.concatenate((levels, levels2))
    levels.sort()

    coords = mesh.coordinates()
    triangles = []
    for cell in df.cells(mesh):
        cell_vertices = cell.entities(0)
        triangles.append(cell_vertices)
    triangulation = Triangulation(coords[:, 0], coords[:, 1], triangles)
    Z = sf.psi.compute_vertex_values()

    ############################################################################
    # Plot streamlines and velocity distributions

    fig = pyplot.figure(figsize=(15, 2.5))
    ax = fig.add_subplot(111)
    ax.plot([-L1, -L1, 0, 0, L2, L2, -L1], [H1, 0, 0, -H2, -H2, H1, H1], 'k')

    # ax.triplot(triangulation, color='#000000', alpha=0.5, lw=0.25)
    ax.tricontour(
        triangulation, Z, levels, linestyles='solid', colors='#0000AA', linewidths=0.5
    )
    ax.plot(x_inflow + u0_inflow, y_inflow, 'm')
    ax.plot(x_inflow - 6 * (y_inflow - 1) * y_inflow, y_inflow, 'k.', ms=4)

    # ax2 = fig.add_subplot(212)
    # ax2.plot(xf, f)
    # ax2.set_xlim(-L1, L2)

    for xi, _idx, _upcrossing in x1_c0:
        ax.plot(xi, x1_ypos, 'rx')
        # ax.text(xi, x1_ypos+0.1, 'X1=%.3f' % x1, backgroundcolor='white')

    ax.set_title(
        'Re = %g, t = %g, xi = %s' % (Re, time, ', '.join('%.3f' % e[0] for e in x1_c0))
    )
    fig.tight_layout()


def get_probe_values(field, x0, x1, N):
    xvec = numpy.linspace(x0[0], x1[0], N)
    yvec = numpy.linspace(x0[1], x1[1], N)
    zvec = numpy.linspace(x0[2], x1[2], N)

    # Get the value at the probe locations
    res = numpy.array([0.0])
    pos = numpy.array([0.0, 0.0, 0.0])
    probe_values = numpy.zeros_like(xvec)
    for i in range(len(probe_values)):
        pos[:] = (xvec[i], yvec[i], zvec[i])
        field.eval(res, pos)
        probe_values[i] = res[0]

    return xvec, yvec, zvec, probe_values


def get_probe_value(field, pos):
    res = numpy.array([0.0])
    pos = numpy.array(pos)
    field.eval(res, pos)
    return res[0]


def get_zero_crossings(x, y):
    N = len(y)
    res = []
    for i in range(1, N - 1):
        yi = y[i]
        yii = y[i + 1]
        if (yi < 0 and yii >= 0) or (yi >= 0 and yii < 0):
            xi, xii = x[i], x[i + 1]
            dy = yii - yi
            fac = yii / dy
            x0 = fac * xi + (1 - fac) * xii
            res.append((x0, i, yi < yii))
    return res


if __name__ == '__main__':
    import warnings

    warnings.filterwarnings('ignore', 'elementwise comparison failed;')

    import dolfin

    dolfin.set_log_level(dolfin.WARNING)

    for h5_file_name in sys.argv[1:]:
        res = load_simulation(h5_file_name)
        plot_backward_facing_step(res)
    pyplot.show()
