"""
Result file statistics
======================

Show statistics on timestep reports logged during a simulation. Give
the name of one or more log files or restart HDF5 files. This small
command line utility will print average values and standard deviations
along with minimum and maximum values of all logged time series

This is a part of the Ocellaris two-phase solver post-processing
software collection 
"""
from __future__ import print_function
import os
from .files import get_result_file_name
from .results import Results


def get_stats(name, ts):
    avg = ts.mean()
    if avg != 0:
        std = ts.std() / avg * 100
    else:
        std = float('NaN')
    return (name, avg, std, ts.min(), ts.max())


def show_logstats(file_name):
    results = Results(file_name)
    maxlen = 0
    stats = []
    for rname, ts in results.reports.items():
        maxlen = max(maxlen, len(rname))
        stats.append(get_stats(rname, ts))
    
    if 'tstime' in results.reports_x:
        stats.append(get_stats('t', results.reports_x['tstime']))
    
    template = '| %%-%ds' % maxlen
    template2 = template + ' |% 11.3g |% 9.2f %% |% 11.3g |% 11.3g |'
    sep = '+' + '-' * (maxlen + 54) + '+'
    
    print(('%%%ds' % len(sep)) % file_name)
    print(sep)
    print(template % 'Name', '|       Mean |    Std.dev |        Min |        Max |')
    print(sep.replace('-', '='))
    i = 0
    L = 3
    for s in sorted(stats):
        if i % L == 0 and i != 0:
            print(sep)
        print(template2 % s)
        i += 1
    if i % L != 0:
        print(sep)
    print()


def main(args):
    # Get report files to read
    file_names = []
    for fn in args[1:]:
        if os.path.isdir(fn):
            fn = get_result_file_name(fn)
        
        if os.path.isfile(fn):
            file_names.append(fn)
        else:
            print('ERROR: not a file %r' % fn)
            exit(1)
    
    if not file_names:
        print(__doc__[1:])
        print()
        print('ERROR: no result files given!')
    
    else:
        print('Result file statistics')
        print('======================\n')
        
        for fn in file_names:
            show_logstats(fn)


if __name__ == '__main__':
    import sys
    main(sys.argv)
