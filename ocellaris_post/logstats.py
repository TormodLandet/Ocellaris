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
import os, sys, argparse
from .files import get_result_file_name
from .results import Results
from . import textplot


# See https://en.wikipedia.org/wiki/ANSI_escape_code
RED = '\033[91m%s\033[0m'    # ANSI escape code Bright Red
YELLOW = '\033[93m%s\033[0m' # ANSI escape code Bright Yellow
BLUE = '\033[94m%s\033[0m' # ANSI escape code Bright Blue


def info(*argv, color='%s'):
    msg = ' '.join(str(w) for w in argv)
    print(color % msg)


def warn(*argv):
    return info(*argv, color=YELLOW)


def error(*argv):
    return info(*argv, color=RED)


def header(*argv):
    return info(*argv, color=BLUE)


def get_stats(name, ts):
    avg = ts.mean()
    if avg != 0:
        std = ts.std() / avg * 100
    else:
        std = float('NaN')
    return (name, avg, std, ts.min(), ts.max())


def print_stats_table(title, stats, maxlen):
    template = '| %%-%ds' % maxlen
    template2 = template + ' |% 11.3g |% 9.2f %% |% 11.3g |% 11.3g |'
    sep = '+' + '-' * (maxlen + 54) + '+'
    
    header(('%%%ds' % len(sep)) % title)
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
    print(sep)
    

def show_logstats(results, skip_first=0):
    maxlen = 1
    stats = []
    for rname, ts in results.reports.items():
        stats.append(get_stats(rname, ts[skip_first:]))
        maxlen = max(maxlen, len(stats[-1][0]))
    
    if 'tstime' in results.reports_x:
        stats.append(get_stats('t', results.reports_x['tstime']))
    
    print_stats_table(results.file_name, stats, maxlen)
    
    if None not in (results.ndofs, results.ncpus):
        print('    Running on %d CPUs with a total of %d DOFs (%d per CPU)'
              % (results.ncpus, results.ncpus, results.ndofs / results.ncpus))


def show_comparison(all_results, ts_name, skip_first=0):
    # Find common file name prefix
    if len(all_results) > 1:
        prefix = os.path.commonpath([res.file_name for res in all_results])
    else:
        prefix, _ = os.path.split(all_results[0].file_name)
    
    maxlen = 1
    stats = []
    for res in all_results:
        if ts_name in res.reports:
            ts = res.reports[ts_name]
        elif ts_name == 't' and 'tstime' in res.reports_x:
            ts = res.reports_x['tstime']
        else:
            warn('\nWARNING: the time series %r is not present %s' %
                 (ts_name, res.file_name))
            continue
        
        # Get a short name with CPU number (for comparing timings)
        name = res.file_name[len(prefix):]
        if name.startswith('/'):
            name = name[1:]
        name = '%s (%r cpus)' % (name, res.ncpus)
        
        stats.append(get_stats(name, ts[skip_first:]))
        maxlen = max(maxlen, len(stats[-1][0]))
    
    print_stats_table('Statistics for time series %r' % ts_name, stats, maxlen)
    print()


def plot_ts(results, ts_names, skip_first=0, title=None):
    if title is not None:
        header(title)
    
    for ts_name in ts_names:
        if ts_name in results.reports:
            x = results.reports_x[ts_name]
            y = results.reports[ts_name]
        elif ts_name == 't' and 'tstime' in results.reports_x:
            x = results.reports_x['tstime']
            y = results.reports_x['tstime']
        else:
            error('\nERROR: the time series %r is not present %s' %
                  (ts_name, results.file_name))
            sys.exit(1)
        header('\nPlot of time series %r:' % ts_name)
        textplot.plot(x[skip_first:], y[skip_first:], figsize=(100, 15))
    print()


def parse_args(argv):
    """
    Parse arguments
    """
    parser = argparse.ArgumentParser(prog='ocellaris_logstats',
                                     description='Ocellaris log statistics post processor')
    parser.add_argument('resfile', nargs='+',
                        help='Name of log or restartfile')
    parser.add_argument('--skip-first', '-s', type=int, metavar='N', default=0,
                        help='Skip the first N data points')
    parser.add_argument('--plot', '-p', metavar='TSNAME', action='append', default=[],
                        help='Plot the given time series (command line output only)')
    parser.add_argument('--restrict', '-r', metavar='TSNAME', action='append', default=[],
                        help='Only show info about the selected time series')
    return parser.parse_args(argv)


def main(args=None):
    # Parse command line arguments
    if args is None:
        args = sys.argv[1:]
    args = parse_args(args)
    
    # Get report files to read
    file_names = []
    for fn in args.resfile:
        if os.path.isdir(fn):
            fn = get_result_file_name(fn)
        
        if os.path.isfile(fn):
            file_names.append(fn)
        else:
            error('\nERROR: not a file %r' % fn)
            exit(1)
    
    if not file_names:
        print(__doc__[1:])
        print()
        error('\nERROR: no result files given!')
        exit(1)
    
    
    print('Ocellaris result file statistics')
    print('================================')

    allres = []
    for fn in file_names:
        # Load result file
        res = Results(fn)
        allres.append(res)
        title = res.file_name
        
        # Show full info
        if not args.restrict:
            show_logstats(res, args.skip_first)
            title = None
        
        # Plot time series
        if args.plot:
            plot_ts(res, args.plot, args.skip_first, title)
        
    for ts_name in args.restrict:
        show_comparison(allres, ts_name, args.skip_first)


if __name__ == '__main__':
    main()
