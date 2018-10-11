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
import sys
import argparse
import glob
import re
import shlex
from subprocess import Popen, PIPE
from time import sleep, time
import signal
from fcntl import fcntl, F_GETFL, F_SETFL
from os import O_NONBLOCK, read, environ
from ocellaris_post import read_yaml_input_file


DESCRIPTION = __doc__


# Restore signals in non-interactive background shells
signal.signal(signal.SIGINT, signal.default_int_handler)
signal.signal(signal.SIGTERM, signal.default_int_handler)
signal.signal(signal.SIGQUIT, signal.default_int_handler)


# Defaults
DEFAULT_TIMEOUT = 60 * 10
DEFAULT_INTERVAL = 10
DEFAULT_MAX_RESTARTS = 2
DEFAULT_KILL_RETRIES = 10
DEFAULT_KILL_WAIT = 10
DEFAULT_SIGINT_WAIT = 60


# ANSI escape sequence to invert foreground and background
INVERTED = '\033[7m%s\033[27m'
RED = '\033[91m%s\033[0m'  # ANSI escape code Bright Red
YELLOW = '\033[93m%s\033[0m'  # ANSI escape code Bright Yellow
BLUE = '\033[94m%s\033[0m'  # ANSI escape code Bright Blue


def say(text):
    sys.stderr.write(INVERTED % text)
    sys.stderr.flush()


def error(text):
    say(RED % text)


def warn(text):
    say(YELLOW % text)


def info(text):
    say(BLUE % text)


def terminate_simulation(
    p, signum=signal.SIGTERM, wait=DEFAULT_KILL_WAIT, retries=DEFAULT_KILL_RETRIES, silent=False
):
    # How long to try the signal before going into TERM / KILL mode
    # Currently hardcoded to minimum 1 and maximum 3 times
    nsig = max(1, min(3, retries / 2))

    # Try to stop the process
    for i in range(retries):
        # This gives the return code if the process has stopped, othewise None
        if p.poll() is not None:
            break

        # First we try to send the signal, after some tries we just term/kill
        if i < nsig:
            info('Sending signal %d to PID %d\n' % (signum, p.pid))
            p.send_signal(signum)
        elif i == nsig:
            info('Terminating PID %d\n' % p.pid)
            p.terminate()
        else:
            warn('Killing PID %d\n' % p.pid)
            p.kill()

        # Wait for process to handle signal
        t = 0
        while t < wait:
            t += 1
            sleep(1)
            read_process_output(p, silent)
            if p.poll() is not None:
                break
    else:
        error('\nERROR: could not terminate simulation!\n\n')
        return False

    info('\nProcess %r exited with status %r\n\n' % (p.pid, p.poll()))
    return True


def read_process_output(p, silent):
    """
    Return a boolean indicating if the process has output on stdout or not
    IMPORTANT: this requires that O_NONBLOCK has been set on p.stdout!
    """
    got_data = False
    while 1:
        try:
            data = read(p.stdout.fileno(), 10000)
        except OSError:
            data = None
        if not data:
            break
        got_data = True
        if not silent:
            data = data.decode('utf8', 'replace')
            sys.stdout.write(data)
            sys.stdout.flush()

    return got_data


def run_simulation(inp_file, ncpus, interval, timeout, silent, mpirun, pystuck):
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
        runner = mpirun
    elif ncpus > 1:
        runner = mpirun + ['-np', str(ncpus)]

    cmd = runner + ['python3', '-m', 'ocellaris', inp_file]
    if pystuck:
        cmd.append('--pystuck')
    p = Popen(cmd, stdout=PIPE)

    # Make sure we can read from p.stdout without blocking
    flags = fcntl(p.stdout, F_GETFL)
    fcntl(p.stdout, F_SETFL, flags | O_NONBLOCK)

    try:
        last_io = time()
        while True:
            exit_code = p.poll()
            if exit_code is not None:
                read_process_output(p, silent)
                return ('exited', exit_code)

            now = time()
            has_output = read_process_output(p, silent)
            if has_output:
                last_io = now
            else:
                # No data to read
                time_since_last_io = now - last_io
                if time_since_last_io > timeout:
                    warn(
                        '\nStdout timeout exceeded, %d seconds since last output\n'
                        % time_since_last_io
                    )

                    if pystuck:
                        info('@@@@@@@@@@ Using pystuck to investigate hung process @@@@@@@@@@')
                        import pystuck

                        pystuck.run_client(stacks=True, ipython=False)
                        info('@@@@@@@@@@ End of pystuck output @@@@@@@@@@')

                    info('Killing child process with PID %d\n' % p.pid)
                    term_ok = terminate_simulation(p, silent=silent)
                    read_process_output(p, silent)
                    return ('timeout', term_ok)
            sleep(interval)
    except KeyboardInterrupt:
        # We got SIGINT, tell Ocellaris to stop
        warn('\nGot SIGINT, letting Ocellaris save restart file\n')
        term_ok = terminate_simulation(
            p, signum=signal.SIGSTOP, wait=DEFAULT_SIGINT_WAIT, silent=silent
        )
        read_process_output(p, silent)
        return ('exited', -2)


def get_restart_files(inp_file):
    """
    Get a list of files that can be used to run the simulation. The first
    returned item is allways the input file name, then follows the names of any
    save point files in order such that the latest version is always last in
    the returned list
    """
    inp = read_yaml_input_file(inp_file)

    prefix = inp['output']['prefix']
    restart_files = glob.glob(prefix + '_savepoint_*.h5')
    restart_files = [
        rf for rf in restart_files if re.match('.*_savepoint_[0-9]+.h5', rf) is not None
    ]

    return [inp_file] + sorted(restart_files)


def babysit_simulation(inp_file, ncpus, interval, timeout, max_restarts, silent, mpirun, pystuck):
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
        info('Running Ocellaris simulation number %d with input file %r\n' % (restarts, fn))

        res, status = run_simulation(fn, ncpus, interval, timeout, silent, mpirun, pystuck)
        if res in 'exited':
            sys.exit(status)
        elif not status:
            error('\nERROR: Giving up!\n')
            sys.exit(31)


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog='orun', description=DESCRIPTION, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input_file', help='Name of inputfile on YAML format')
    parser.add_argument(
        '--ncpus', '-n', metavar='NCPU', default=1, type=int, help='Number of MPI processes'
    )
    parser.add_argument(
        '--interval',
        '-i',
        metavar='INTERVAL',
        type=float,
        default=DEFAULT_INTERVAL,
        help='Output interval in seconds',
    )
    parser.add_argument(
        '--pystuck', action='store_true', help='Enable pystuck on the root MPI rank.'
    )
    parser.add_argument(
        '--timeout',
        '-t',
        metavar='TIMEOUT',
        type=float,
        default=DEFAULT_TIMEOUT,
        help='Output timeout in seconds. After this period of inactivity '
        'the simulation is killed',
    )
    parser.add_argument(
        '--restarts',
        '-r',
        metavar='RESTARTS',
        type=int,
        default=DEFAULT_MAX_RESTARTS,
        help='Number of restarts of the same file (input or restart file). '
        'Every time the simulation writes a new savepoint the counter is reset',
    )
    parser.add_argument(
        '--silent',
        '-s',
        action='store_true',
        default=False,
        help='Do not relay stdout from Ocellaris',
    )
    parser.add_argument(
        '--mpirun', metavar='MPIRUN', default='mpirun', help='The mpirun executable'
    )
    return parser.parse_args(argv)


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    args = parse_args(args)
    babysit_simulation(
        inp_file=args.input_file,
        ncpus=args.ncpus,
        interval=args.interval,
        timeout=args.timeout,
        max_restarts=args.restarts,
        silent=args.silent,
        mpirun=shlex.split(args.mpirun),
        pystuck=args.pystuck,
    )


if __name__ == '__main__':
    main()
