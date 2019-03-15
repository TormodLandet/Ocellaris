# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
import numpy


def get_root_value(value, comm=None):
    """
    Return the value that is given on the root process
    """
    ncpu = dolfin.MPI.size(dolfin.MPI.comm_world)

    if ncpu == 1:
        # Not running in parallel
        return value

    if comm is None:
        comm = dolfin.MPI.comm_world  # a mpi4py communicator
    return comm.bcast(value)


def sync_arrays(array_list, sync_all=False, comm=None):
    """
    Given a list of arrays on each process (dtype=float, various length
    and number of arrays on each process), make sure all arrays end up
    on the root process

    If sync_all is true then all processes end up with the same data
    """
    if comm is None:
        comm = dolfin.MPI.comm_world  # a mpi4py communicator

    # Check if we are running in parallel or not
    rank = dolfin.MPI.rank(comm)
    ncpu = dolfin.MPI.size(comm)
    if ncpu == 1:
        # All arrays are allready on the root process
        return

    # Receive on root (rank 0), send on all other ranks > 0
    if rank == 0:
        # Loop through non-root processes and get their lines
        for proc in range(1, ncpu):
            _sync_array_list_one_way(
                comm, rank, array_list, from_rank=proc, to_rank=rank
            )
    else:
        _sync_array_list_one_way(comm, rank, array_list, from_rank=rank, to_rank=0)

    # Syncronize all processes
    comm.barrier()

    if not sync_all:
        return

    # Send on root rank and recieve on all other ranks
    if rank == 0:
        # Loop through non-root processes and send the result to them
        for proc in range(1, ncpu):
            _sync_array_list_one_way(
                comm, rank, array_list, from_rank=rank, to_rank=proc
            )
    else:
        del array_list[:]
        _sync_array_list_one_way(comm, rank, array_list, from_rank=0, to_rank=rank)

    # Syncronize all processes
    comm.barrier()


def _sync_array_list_one_way(comm, rank, array_list, from_rank, to_rank):
    """
    Send a list of arrays from one rank to be appended in the list
    of arrays in a different rank
    """
    if rank == from_rank:
        # Send number of arrays
        N = len(array_list)
        comm.isend(N, dest=to_rank)

        if N > 0:
            # Send array lengths
            array_lengths = numpy.empty(N, int)
            for i, arr in enumerate(array_list):
                M, = arr.shape  # implicitly checks shape
                array_lengths[i] = M
            comm.Isend(array_lengths, dest=to_rank)

            # Send arrays
            for arr in array_list:
                comm.Isend(arr, dest=to_rank)

    else:
        assert rank == to_rank

        # Recieve number of arrays
        N = comm.irecv(source=from_rank).wait()

        if N > 0:
            # Recieve line lengths
            array_lengths = numpy.empty(N, int)
            comm.Irecv(array_lengths, source=from_rank).Wait()

            # Recieve arrays
            for M in array_lengths:
                arr = numpy.empty(M, float)
                comm.Irecv(arr, source=from_rank).Wait()
                array_list.append(arr)


def gather_lines_on_root(lines, comm=None):
    """
    Given a list of lines, use MPI to add all processes' lines to the list of
    lines on the root process.

    The list of lines MUST be a list of (x, y) tuples where x and y are numpy
    arrays with dtype=float and len(x) == len(y)
    """
    if comm is None:
        comm = dolfin.MPI.comm_world  # a mpi4py communicator

    # Check if we are running in parallel or not
    if dolfin.MPI.size(comm) == 1:
        # All arrays are allready on the root process
        return

    # Split point arrays coordinates into separate items
    array_list = []
    for x, y in lines:
        array_list.append(x)
        array_list.append(y)

    # Get arrays onto root process
    sync_arrays(array_list, sync_all=False, comm=comm)

    # Merge point arrays coordinates
    del lines[:]
    N = len(array_list)
    assert N % 2 == 0
    for i in range(0, N, 2):
        lines.append((array_list[i], array_list[i + 1]))
