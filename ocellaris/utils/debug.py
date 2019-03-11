# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import sys
import dolfin


FUNC_TRACE_PATTERNS = ['ocellaris', 'solenoidal']
LINE_TRACE_PATTERNS = []


def enable_super_debug(stderr=False, func_trace_patterns=None, line_trace_patterns=None):
    """
    For those times when a C++ call is hanging and you need to know where,
    or the code SEGFAULTs without any sensible output. This function installs
    a execution tracer into the Python runtime which logs all function calls
    and possibly some line for line running logs for files you really want to
    interrogate.
    
    This will produce a LOT of output and slow down the program significantly,
    but sometimes it is the easiest option. When running in pytest try something
    like::

        mpirun -np 4 python3 -m pytest -s  --instafail ...

    This will show the exception instead of just hanging on the MPI barriers
    defined in tests/conftest.py

    NEVER include the Python standard library in ``func_trace_patterns`` or
    ``line_trace_patterns``, otherwise the act of logging generate several
    function calls to be logged and so on which tends to lock up rather hard.
    These lists contain (parts of) full ``*.py`` file paths that you want to
    investigate by this excessive debugg logging functionality.
    """
    # What should we trace
    if func_trace_patterns is None:
        func_trace_patterns = FUNC_TRACE_PATTERNS
    if line_trace_patterns is None:
        line_trace_patterns = line_trace_patterns

    rank = dolfin.MPI.rank(dolfin.MPI.comm_world)

    # Open output file for writing the debug log
    if stderr:
        out_name = 'sys.stderr'
        outfile = sys.stderr
    else:
        out_name = 'OCELLARISSUPERDEBUG_%d' % rank
        outfile = open(out_name, 'wt')

    def trace_lines(frame, event, arg):
        """
        Print every line that is about to run
        """
        if event != 'line':
            return
        co = frame.f_code
        func_name = co.co_name
        line_no = frame.f_lineno
        outfile.write('   Rank %d about to run %s line %s\n' % (rank, func_name, line_no))
        outfile.flush()

    def trace(frame, event, arg):
        """
        Print every function/method call that is about to happen
        """
        if event != 'call':
            return

        func_name = frame.f_code.co_name
        file_name = frame.f_code.co_filename
        line_no = frame.f_lineno

        # Ignore deep dives in other libraries
        caller = frame.f_back
        if caller is None:
            return
        caller_name = caller.f_code.co_name
        caller_file_name = caller.f_code.co_filename
        caller_line_no = caller.f_lineno
        want_to_know = False
        for interesting in func_trace_patterns:
            if interesting in caller_file_name:
                want_to_know = True

        if want_to_know:
            # Must NEVER run when stdout.write() or flush() will trigger this trace
            outfile.write(
                'Rank %d has call to %s (%s @ %s) from %s (%s @ %s)\n'
                % (
                    rank,
                    func_name,
                    file_name,
                    line_no,
                    caller_name,
                    caller_file_name,
                    caller_line_no,
                )
            )
            outfile.flush()

            # Trace the lines of the function that is about to run?
            for interesting in line_trace_patterns:
                if interesting in file_name:
                    return trace_lines

    sys.stdout.write('Enabling SUPER DEBUG - logging to %s\n' % out_name)
    sys.settrace(trace)
