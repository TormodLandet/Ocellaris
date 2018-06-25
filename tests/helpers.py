import os
import dolfin
import pytest


# -- MPI helpers ---------------------------------------------------------


COMM = dolfin.MPI.comm_world
RANK = COMM.rank
SIZE = COMM.size
skip_in_parallel = pytest.mark.skipif(
    SIZE > 1, reason='This test should ' 'only be run in serial.'
)


def all_ok(status):
    """
    Given a boolean status on each rank, return whether the status is OK
    on all MPI ranks or not (boolean)
    """
    stat = 1.0 if status else 0.0
    all_stat = dolfin.MPI.min(dolfin.MPI.comm_world, stat)
    return all_stat > 0


def mpi_print(arg1, *args, **kwargs):
    if arg1.startswith('\n'):
        arg1 = '\nRANK%d/%d: %s' % (RANK, SIZE, arg1[1:])
    else:
        arg1 = 'RANK%d/%d: %s' % (RANK, SIZE, arg1)
    print(arg1, *args, **kwargs)


def mpi_int_sum(value):
    v = dolfin.MPI.sum(dolfin.MPI.comm_world, float(value))
    return int(round(v))


def comm_self_to_comm_world(s, w):
    """
    Convert a function defined on a mesh with a MPI_COMM_SELF to
    a function defined on the same geometry, but with MPI_COMM_WORLD.
    The mesh and function spaces mush be identical except for the
    MPI communicator.

    - s is the MPI_COMM_SELF function
    - w is the MPI_COMM_WORLD function

    The return value is in the same function space as w
    """
    Vs = s.function_space()
    Vw = w.function_space()
    gdim = Vs.mesh().geometry().dim()
    coords_s = Vs.tabulate_dof_coordinates().reshape((-1, gdim))
    coords_w = Vw.tabulate_dof_coordinates().reshape((-1, gdim))

    def locate_in_s(coord):
        for i, c in enumerate(coords_s):
            if (c == coord).all():
                return i

    # Construct w2 and initialize to zero
    w2 = w.copy()
    w2.vector().zero()
    w2.vector().apply('insert')

    # Construct a parallel version of s by finding dofs with the
    # exact same geometrical coordinates
    arr_s = s.vector().get_local()
    arr_w2 = w2.vector().get_local()
    for dof_w, coord in enumerate(coords_w):
        dof_s = locate_in_s(coord)
        arr_w2[dof_w] = arr_s[dof_s]
    w2.vector().set_local(arr_w2)
    w2.vector().apply('insert')
    return w2


# -- File handling -------------------------------------------------------


def get_test_file_name(fn):
    mydir = os.path.dirname(__file__)
    return os.path.join(mydir, 'data', fn)


def mpi_tmpdir(tmpdir_factory, dir_name):
    """
    Returns the same directory name on all MPI processes
    """
    comm = dolfin.MPI.comm_world
    if comm.rank == 0:
        x = tmpdir_factory.mktemp(dir_name, numbered=True)
        x = str(x)
    else:
        x = None
    return comm.bcast(x)
