"""
Read a dolfin h5 file (Ocellaris restart file) and extract a 2D slice
of velocities and pressures which is stored as a numpy array on disk

Assumes that the simulation is 2D (ignores the z direction)
"""
import numpy
import dolfin as df


Nx = Ny = 20
startx, endx = 0, 5.5
starty, endy = -0.9423, 1


def load_data(h5_file_name):
    """
    Use the dolfin HDF5 classes to load velocity and
    pressure fields from the h5 file along with the
    mesh
    """
    h5 = df.HDF5File(df.MPI.comm_world, h5_file_name, 'r')

    mesh = df.Mesh()
    h5.read(mesh, '/mesh', False)

    # Read velocities
    u0_signature = h5.attributes('/u0')['signature']
    eu = eval('df.' + u0_signature.replace('triangle', 'df.triangle'))
    Vu = df.FunctionSpace(mesh, eu)
    u0 = df.Function(Vu)
    u1 = df.Function(Vu)
    h5.read(u0, '/u0')
    h5.read(u1, '/u1')

    # Read pressure
    p_signature = h5.attributes('/p')['signature']
    ep = eval('df.' + p_signature.replace('triangle', 'df.triangle'))
    Vp = df.FunctionSpace(mesh, ep)
    p = df.Function(Vp)
    h5.read(p, '/p')

    h5.close()
    res = {'u0': u0, 'u1': u1, 'p': p}

    return res


def get_field_slice(field, startx, endx, starty, endy, Nx, Ny):
    """
    Extract a regular spaced slice of the field. Returns
    full coordinate and value matrices (numpy meshgrid style)
    """
    xvec = numpy.linspace(startx, endx, Nx)
    yvec = numpy.linspace(starty, endy, Ny)

    X, Y = numpy.meshgrid(xvec, yvec)
    V = numpy.zeros_like(X)
    res = numpy.array([0.0])
    pos = numpy.array([0.0, 0.0, 0.0])
    for i in xrange(Nx):
        for j in xrange(Ny):
            pos[:2] = (X[i, j], Y[i, j])
            field.eval(res, pos)
            V[i, j] = res[0]

    return X, Y, V


def save_ndarrays(h5_file_name, numpy_file_name, startx, endx, starty, endy, Nx, Ny):
    """
    Extract velocity and pressure fields and save to a *.npy file
    """
    data = load_data(h5_file_name)
    X, Y, U0 = get_field_slice(data['u0'], startx, endx, starty, endy, Nx, Ny)
    X, Y, U1 = get_field_slice(data['u1'], startx, endx, starty, endy, Nx, Ny)
    X, Y, P = get_field_slice(data['p'], startx, endx, starty, endy, Nx, Ny)

    numpy.save(numpy_file_name, numpy.array([X, Y, U0, U1, P]))


if __name__ == '__main__':
    import sys

    h5_file_name = sys.argv[1]
    numpy_file_name = sys.argv[2]

    # Load h5 data and save to numpy data format
    save_ndarrays(h5_file_name, numpy_file_name, startx, endx, starty, endy, Nx, Ny)

    # Load numpy data format and plot the loaded fields
    X, Y, U, V, P = numpy.load(numpy_file_name)
    from matplotlib import pyplot

    pyplot.contourf(X, Y, P)
    pyplot.quiver(X, Y, U, V)
    pyplot.show()
