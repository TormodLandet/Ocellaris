import dolfin
import numpy


def get_root_value(value):
    """
    Return the value that is given on the root process
    """
    ncpu = dolfin.MPI.size(dolfin.mpi_comm_world())
    
    if ncpu == 1:
        # Not running in parallel
        return value
    
    from mpi4py.MPI import COMM_WORLD as comm
    return comm.bcast(value)


def gather_lines_on_root(lines):
    """
    Given a list of lines, use MPI to add all processes' lines to the list of
    lines on the root process.
    
    The list of lines MUST be a list of (x, y) tuples where x and y are numpy
    arrays with dtype=float and len(x) == len(y)
    """
    # Check if we are running in parallel or not
    rank = dolfin.MPI.rank(dolfin.mpi_comm_world())
    ncpu = dolfin.MPI.size(dolfin.mpi_comm_world())
    if ncpu == 1:
        # All lines are allready on the root process
        return
    
    from mpi4py.MPI import COMM_WORLD as comm

    # Receive on root (rank 0), send on ranks > 0
    if rank == 0:
        # Loop through non-root processes and get their lines
        for proc in range(1, ncpu):
            # Get number of lines
            num_proc_lines = comm.irecv(source=proc).wait()
            
            if num_proc_lines == 0:
                continue
            
            # Get line lengths
            proc_line_lengths = numpy.empty(num_proc_lines, int)
            comm.Irecv(proc_line_lengths, source=proc).Wait()
            
            for n in proc_line_lengths:
                x = numpy.empty(n, float)
                y = numpy.empty(n, float)
                comm.Irecv(x, source=proc).Wait()
                comm.Irecv(y, source=proc).Wait()
                lines.append((x, y))
        
    else:
        # Send number of lines to the root process
        comm.isend(len(lines), dest=0)
        
        if len(lines) > 0:
            # Send line lengths
            line_lengths = numpy.array([len(line[0]) for line in lines], int)
            comm.Isend(line_lengths, dest=0)
            
            # Send lines
            for x, y in lines:
                comm.Isend(x, dest=0)
                comm.Isend(y, dest=0)
    
    # Syncronize all processes
    comm.barrier()


def test_gather_lines_on_root():
    from matplotlib import pyplot
    import random
    
    # Create a random number of random lines
    num_lines = random.randint(0, 6)    
    lines = []
    for _ in range(num_lines):
        len_line = random.randint(1, 10)
        x = numpy.array([random.random() for _ in range(len_line)], float)
        y = numpy.array([random.random() for _ in range(len_line)], float)
        lines.append((x, y))

    # Send lines to rank 0
    gather_lines_on_root(lines)

    # Plot the process lines (on ranks > 0) or all lines (on rank 0)
    pyplot.figure()
    for x, y in lines:
        pyplot.plot(x, y)
    pyplot.xlim(0, 1)
    pyplot.ylim(0, 1)
    pyplot.show()


if __name__ == '__main__':
    test_gather_lines_on_root()
