import numpy
import dolfin
from ocellaris.utils import gather_lines_on_root, sync_arrays


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
