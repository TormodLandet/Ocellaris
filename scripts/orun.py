#!/usr/bin/env python3
"""
orun
====

Start an Ocellaris simulation with a "babysitter" that watches the stdout
stream and kills the simulation if no output is produced over an extended
period of time (default 10 minutes / 600 seconds).

If the simulation writes restart files at regular intervals then the babysitter
can be made to restart the simulation from the latest restart file a number of
times (default 2 restarts of the same file).

The reason for this babysitter is that there are difficult to debug problems
(probably in PETSc) that causes the simulation to be stuck at 100% CPU
utilisation. No backtrace is available on Ctrl+C / SIGINT which would be the
case if there was a infinite loop in the Python code, so most likely the error
exists in a C extension.
"""
DESCRIPTION = __doc__

import sys
import argparse
import glob
import re
from subprocess import Popen, PIPE
from time import sleep, time
import signal
from fcntl import fcntl, F_GETFL, F_SETFL
from os import O_NONBLOCK, read, environ
import yaml


# Restore signals in non-interactive background shells
signal.signal(signal.SIGINT, signal.default_int_handler)
signal.signal(signal.SIGTERM, signal.default_int_handler)
signal.signal(signal.SIGQUIT, signal.default_int_handler)


# Defaults
DEFAULT_TIMEOUT = 60 * 10
DEFAULT_INTERVAL = 10
DEFAULT_MAX_RESTARTS = 2
DEFAULT_KILL_TRIES = 10
DEFAULT_KILL_WAIT = 10
DEFAULT_SIGINT_WAIT = 60


# ANSI escape sequence to invert foreground and background
INVERTED = '\033[7m%s\033[27m'
RED = '\033[91m%s\033[0m'    # ANSI escape code Bright Red
YELLOW = '\033[93m%s\033[0m' # ANSI escape code Bright Yellow
BLUE = '\033[94m%s\033[0m' # ANSI escape code Bright Blue


def say(text):
    sys.stderr.write(INVERTED % text)
    sys.stderr.flush()


def error(text):
    say(RED % text)


def warn(text):
    say(YELLOW % text)


def info(text):
    say(BLUE % text)


def terminate_simulation(p, signal=signal.SIGTERM,
                         wait=DEFAULT_KILL_WAIT,
                         retries=DEFAULT_KILL_TRIES):
    p.send_signal(signal)
    for _ in range(retries):
        sleep(wait)
        if p.poll() is not None:
            break
        p.kill()
    else:
        error('\nERROR: could not terminate simulation!\n\n')
        return False
    return True


def run_simulation(inp_file, ncpus, interval, timeout, silent):
    """
    Run the simulation while watching the stdout for output every ``interval``
    seconds and terminating after ``timeout`` seconds without any output.
    
    Returns one of ``('exited', exit_code)``  or ``('timeout', term_ok)``.
    
    Where ``term_ok`` is ``True`` if the simulation was shut down properly
    after the timeout and ``False`` if it is still running (unkillable). This
    should probably never happen, but you never know with IO drivers and
    other comblex networked systems... This way we avoid running more than one
    instance of the simulation at any time.
    """
    runner = []
    in_queue = environ.get('SLURM_NTASKS', False)
    if in_queue:
        runner = ['mpirun']
    elif ncpus > 1:
        runner = ['mpirun', '-np',  str(ncpus)]
    
    cmd = runner + ['python3', '-m', 'ocellaris', inp_file]
    p = Popen(cmd, stdout=PIPE)

    # Make sure we can read from p.stdout without blocking
    flags = fcntl(p.stdout, F_GETFL)
    fcntl(p.stdout, F_SETFL, flags | O_NONBLOCK)
    
    try:
        last_io = time()
        while True:
            exit_code = p.poll()
            if exit_code is not None:
                return ('exited', exit_code)
            
            now = time()
            try:
                data = read(p.stdout.fileno(), 10000)
                if not silent:
                    data = data.decode('utf8', 'replace')
                    sys.stdout.write(data)
                    sys.stdout.flush()
                last_io = now
            except OSError:
                # No data to read
                time_since_last_io = now - last_io 
                if time_since_last_io > timeout:
                    warn('\nStdout timeout exceeded, %d seconds since last output\n'
                         % time_since_last_io)
                    info('Killing child process with PID %d\n' % p.pid)
                    term_ok = terminate_simulation(p)
                    return ('timeout', term_ok)
            sleep(interval)
    except KeyboardInterrupt:
        # We got SIGINT, tell Ocellaris to stop
        warn('\nGot SIGINT, letting Ocellaris save restart file\n')
        info('Stopping child process with PID %d\n' % p.pid)
        term_ok = terminate_simulation(p, signal=signal.SIGINT, wait=DEFAULT_SIGINT_WAIT)
        return ('exited', -2)


def get_restart_files(inp_file):
    """
    Get a list of files that can be used to run the simulation. The first
    returned item is allways the input file name, then follows the names of any
    save point files in order such that the latest version is always last in
    the returned list
    """
    with open(inp_file, 'rt') as inp:
        input_str = inp.read()
    inp = yaml.load(input_str)
    
    prefix = inp['output']['prefix']
    restart_files = glob.glob(prefix + '_savepoint_*.h5')
    restart_files = [rf for rf in restart_files if re.match('.*_savepoint_[0-9]+.h5', rf) is not None]
    
    return [inp_file] + sorted(restart_files)


def babysit_simulation(inp_file, ncpus, interval, timeout, max_restarts, silent):
    """
    Run the simulation, possibly restarting from savepoint files in case of a
    stuck simulation. Only ``max_restarts`` of the same file are tried before
    giving up and exiting with an error
    """
    restarts = 0
    nfiles = 0
    while True:
        files = get_restart_files(inp_file)
        
        # Manage restart limit
        if len(files) > nfiles:
            restarts = 0
            nfiles = len(files)
        restarts += 1
        
        if restarts >= max_restarts + 2:
            error('\nMaximum number of restarts exceeded\n')
            return
         
        fn = files[-1]
        info('Running Ocellaris simulation number %d with input file %r\n'
             % (restarts, fn))
        
        res, status = run_simulation(fn, ncpus, interval, timeout, silent)
        if res in 'exited':
            sys.exit(status)
        elif not status:
            error('\nERROR: Giving up!\n')
            sys.exit(31)


def parse_args(argv):
    parser = argparse.ArgumentParser(prog='orun',
                                     description=DESCRIPTION)
    parser.add_argument('input_file', help='Name of inputfile on YAML format')
    parser.add_argument('--ncpus', '-n', metavar='NCPU', default=1, type=int,
                        help='Number of MPI processes')
    parser.add_argument('--interval', '-i', metavar='INTERVAL', type=float, default=DEFAULT_INTERVAL,
                        help='Output interval in seconds')
    parser.add_argument('--timeout', '-t', metavar='TIMEOUT', type=float, default=DEFAULT_TIMEOUT,
                        help='Output timeout in seconds. After this period of inactivity '
                        'the simulation is killed')
    parser.add_argument('--restarts', '-r', metavar='RESTARTS', type=int, default=DEFAULT_MAX_RESTARTS,
                        help='Number of restarts of the same file (input or restart file). '
                        'Every time the simulation writes a new savepoint the counter is reset')
    parser.add_argument('--silent', '-s', action='store_true', default=False,
                        help='Do not relay stdout from Ocellaris')
    return parser.parse_args(argv)


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    args = parse_args(args)
    babysit_simulation(inp_file=args.input_file,
                       ncpus=args.ncpus,
                       interval=args.interval,
                       timeout=args.timeout,
                       max_restarts=args.restarts,
                       silent=args.silent)


if __name__ == '__main__':
    main()
