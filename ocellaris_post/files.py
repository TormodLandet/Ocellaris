import os

def get_result_file_name(result_directory):
    """
    Given a directory of Ocellaris results from one simulation,
    return the best result file name
    """
    assert os.path.isdir(result_directory)
    
    logfiles = []
    endpoints = []
    for fn in os.listdir(result_directory):
        if fn.endswith('.log'):
            logfiles.append(os.path.join(result_directory, fn))
        elif 'endpoint' in fn and fn.endswith('.h5'):
            endpoints.append(os.path.join(result_directory, fn))
    
    if len(endpoints) == 1 and len(logfiles) == 0:
        return endpoints[0]
    
    elif len(endpoints) == 1 and len(logfiles) == 1:
        t1 = os.path.getmtime(endpoints[0])
        t2 = os.path.getmtime(logfiles[0])
        if t2 - t1 > 60:
            raise RuntimeError('Found logfile and endpoint file in %r' % result_directory
                           + ' logfile is significantly newer than endpointfile!')
        return endpoints[0]
    
    elif len(endpoints) > 1 or len(logfiles) > 1:
        raise RuntimeError('Multiple logfiles and or endpoints found in %r' % result_directory
                           + ' could not find unique result file')
    
    elif len(logfiles) == 1:
        return logfiles[0]
    
    else:
        raise RuntimeError('No logfiles and or endpoints found in %r' % result_directory)