from functools import wraps
from collections import defaultdict
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
        timer = dolfin.Timer(task)
        ret =  f(*args, **kwds)
        t = timer.stop()
        timeit.timings[task].append(t)
        return ret
    return wrapper
timeit.timings = defaultdict(list)
