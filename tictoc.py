import time, threading

def tic(description=''):
    data = threading.local()
    if not hasattr(data, 'tictoc'):
        data.tictoc = []
    data.tictoc.append((description, time.time()))
    print 'START "%s"' % description

def toc():
    data = threading.local()
    if not hasattr(data, 'tictoc') or not data.tictoc:
        print 'ERROR: toc() without a tic()'
        return
    
    description, t0 = data.pop() 
    dt = time.time() - t0
    print 'DONE "%s" in %.2fs' % (description, dt)
