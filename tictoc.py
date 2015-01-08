import time

_tictoc_data = []

def tic(description=''):
    _tictoc_data.append((description, time.time()))
    print 'START "%s"' % description

def toc():
    if not _tictoc_data:
        print 'ERROR: toc() without a tic()'
        return
    
    description, t0 = _tictoc_data.pop() 
    dt = time.time() - t0
    print 'DONE "%s" in %.2fs' % (description, dt)

if __name__ == '__main__':
    tic()
    time.sleep(0.5)
    toc()
    
    tic('test2')
    toc()
