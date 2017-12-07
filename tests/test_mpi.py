import numpy
import dolfin
from ocellaris.utils import gather_lines_on_root, sync_arrays, get_local, set_local


def test_gather_points_on_root():
    comm = dolfin.MPI.comm_world # a mpi4py communicator
    rank = dolfin.MPI.rank(comm)
    ncpu = dolfin.MPI.size(comm)
    num_lines = 3
    len_line = 10
    
    lines = []
    for _ in range(num_lines):
        x = numpy.arange(len_line, dtype=float) * rank
        y = numpy.arange(len_line, dtype=float) * rank * 2
        lines.append((x, y))
    
    # Send lines to rank 0
    gather_lines_on_root(lines)
    
    def check_lines(lines, r):
        assert len(lines) == num_lines
        for x, y in lines:
            assert x.shape == y.shape == (len_line,)
            assert (x == numpy.arange(len_line, dtype=float) * r).all()
            assert (y == numpy.arange(len_line, dtype=float) * r * 2).all()
    
    if rank != 0:
        check_lines(lines, rank)
    else:
        assert len(lines) == num_lines * ncpu
        for i in range(ncpu):
            ls = lines[i*num_lines : (i+1)*num_lines]
            check_lines(ls, i)


def test_sync_arrays():
    comm = dolfin.MPI.comm_world # a mpi4py communicator
    rank = dolfin.MPI.rank(comm)
    ncpu = dolfin.MPI.size(comm)
    
    # Construct some points
    num_points = 20
    dim = 3
    points = []
    for _ in range(num_points):
        points.append(numpy.arange(dim, dtype=float) * rank)
    
    # Sync all points
    sync_arrays(points, sync_all=True, comm=comm)
    
    def check_lines(pts, r):
        assert len(pts) == num_points
        for pt in pts:
            assert pt.shape == (dim,)
            assert (pt == numpy.arange(dim, dtype=float) * r).all()
    
    # Verify that all points ended up on all processes
    assert len(points) == num_points * ncpu
    for i in range(ncpu):
        ls = points[i*num_points : (i+1)*num_points]
        check_lines(ls, i)


def test_get_set_local_with_ghosts():
    comm = dolfin.MPI.comm_world
    rank = dolfin.MPI.rank(comm)
    
    dolfin.parameters['ghost_mode'] = 'shared_vertex'
    mesh = dolfin.UnitSquareMesh(comm, 4, 4)
    V = dolfin.FunctionSpace(mesh, 'DG', 0)
    u = dolfin.Function(V)
    dm = V.dofmap()
    im = dm.index_map()
    
    # Write global dof number into array for all dofs
    start, end = u.vector().local_range()
    arr = u.vector().get_local()
    arr[:] = numpy.arange(start, end)
    u.vector().set_local(arr)
    u.vector().apply('insert')
    
    # Get local + ghost values
    arr2 = get_local(u.vector(), V)
    
    # Get ghost global indices and some sizes
    global_ghost_dofs = im.local_to_global_unowned()
    Nall = im.size(im.MapSize.ALL)
    Nghost = global_ghost_dofs.size
    Nown = im.size(im.MapSize.OWNED)
    assert Nown + Nghost == Nall
    
    # Produce the expected answer
    dofs = numpy.arange(start, end + Nghost)
    dofs[Nown:] = global_ghost_dofs
    
    # Check the results
    assert dofs.shape == arr2.shape
    diff = abs(dofs - arr2).max()
    print
    print(rank, start, end, global_ghost_dofs)
    print(rank, numpy.array(arr2, dtype=numpy.intc), '\n ', dofs, diff)
    assert diff == 0
