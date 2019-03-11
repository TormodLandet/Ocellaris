# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from functools import wraps
import time
import dolfin


def timeit(f):
    """
    Timer decorator

    This decorator stores the cummulative time spent in each
    function that is wrapped by the decorator.

    Functions are identified by their names
    """
    # Extract name at defition time, not at run time
    if timeit.next_name:
        # Used by timeit_named as a way to pass information
        timed_task_name = timeit.next_name
        timeit.next_name = None
    else:
        timed_task_name = f.__name__

    @wraps(f)
    def wrapper(*args, **kwds):
        with dolfin.Timer('Ocellaris %s' % timed_task_name):
            # print('<%s>' % timed_task_name)
            ret = f(*args, **kwds)
            # print('</%s>' % timed_task_name)
        return ret

    return wrapper


timeit.next_name = None


def timeit_named(name):
    """
    Generate a decorator that uses the given name instead on the
    callable functions __name__ in the timings table
    """
    timeit.next_name = name
    return timeit


timeit.named = timeit_named


def log_timings(simulation, clear=False):
    """
    Print the FEniCS + Ocellaris timings to the log
    """
    # Total time spent in the simulation
    tottime = time.time() - simulation.t_start

    # Get timings from FEniCS and sort by total time spent
    tclear = dolfin.TimingClear.clear if clear else dolfin.TimingClear.keep
    timingtypes = [
        dolfin.TimingType.user,
        dolfin.TimingType.system,
        dolfin.TimingType.wall,
    ]
    table = dolfin.timings(tclear, timingtypes)
    table_lines = table.str(True).split('\n')
    simulation.log.info('\nFEniCS timings:   %s  wall pst' % table_lines[0][18:])
    simulation.log.info(table_lines[1] + '-' * 10)
    tmp = [(float(line.split()[-5]), line) for line in table_lines[2:]]
    tmp.sort(reverse=True)
    for wctime, line in tmp:
        simulation.log.info('%s    %5.1f%%' % (line, wctime / tottime * 100))
