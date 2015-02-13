from functools import wraps
import time
from collections import defaultdict

def timeit(f):
    """
    Timer decorator
    
    This decorator stores the cummulative time spent in each 
    function that is wrapped by the decorator.
    
    Functions are identified by their names
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        t1 = time.time()
        ret =  f(*args, **kwds)
        timeit.timings[f.__name__].append(time.time()-t1)
        return ret
    return wrapper
timeit.timings = defaultdict(list)
