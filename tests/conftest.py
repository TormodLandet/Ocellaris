import dolfin


# Get the MPI_COMM_WORLD rank and size
COMM = dolfin.MPI.comm_world
RANK = COMM.rank
SIZE = COMM.size


def pytest_runtest_logstart(nodeid, location):
    """
    Help pinpoint MPI deadlocks by putting a barrier and an stdin flush
    before and after each test run

    This runs before each test
    """
    if SIZE == 1:
        return
    print()
    print('Going to run %s on %d/%d' % (nodeid, RANK, SIZE), flush=True)
    COMM.barrier()
    print('BARRIER_1 done %d/%d' % (RANK, SIZE), flush=True)


def pytest_runtest_logfinish(nodeid, location):
    """
    Help pinpoint MPI deadlocks by putting a barrier and an stdin flush
    before and after each test run

    This runs after each
    """
    if SIZE == 1:
        return
    print()
    print('BARRIER_2 running now on %d/%d' % (RANK, SIZE), flush=True)
    COMM.barrier()
    print('Finished running %s on %d/%d' % (nodeid, RANK, SIZE), flush=True)
