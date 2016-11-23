from functools import wraps
from collections import defaultdict
import time
import dolfin

def timeit(f):
    """
    Timer decorator
    
    This decorator stores the cummulative time spent in each 
    function that is wrapped by the decorator.
    
    Functions are identified by their names
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        task = f.__name__
        timer = dolfin.Timer('Ocellaris %s' % task)
        #print '<%s>' % task
        ret =  f(*args, **kwds)
        #print '</%s>' % task 
        t = timer.stop()
        timeit.timings[task].append(t)
        return ret
    return wrapper
timeit.timings = defaultdict(list)


def log_timings(simulation):
    """
    Print the FEniCS + Ocellaris timings to the log
    """
    # Total time spent in the simulation
    tottime = time.time() - simulation.t_start
    
    # Get timings from FEniCS and sort by total time spent
    timingtypes = [dolfin.TimingType_user, dolfin.TimingType_system, dolfin.TimingType_wall]
    table = dolfin.timings(dolfin.TimingClear_keep, timingtypes)
    table_lines = table.str(True).split('\n')
    simulation.log.info('\nFEniCS timings:   %s  wall pst' % table_lines[0][18:])
    simulation.log.info(table_lines[1] + '-'*10)
    tmp = [(float(line.split()[-5]), line) for line in table_lines[2:]]
    tmp.sort(reverse=True)
    for wctime, line in tmp:
        simulation.log.info('%s     %4.1f%%' % (line, wctime/tottime*100))

