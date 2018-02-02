import sys
import dolfin


FUNC_TRACE_PATTERNS = ['ocellaris', 'solenoidal']
LINE_TRACE_PATTERNS = []


def enable_super_debug():
    """
    For those times when a C++ call is hanging and you need to know where ...
    Logs all function calls to stderr, will produce a LOT of output and
    slow down the program significantly
    """
    rank = dolfin.MPI.rank(dolfin.MPI.comm_world)
    outfile = open('OCELLARISSUPERDEBUG_%d' % rank, 'wt')
    
    def trace_lines(frame, event, arg):
        if event != 'line':
            return
        co = frame.f_code
        func_name = co.co_name
        line_no = frame.f_lineno
        outfile.write('   About to run %s line %s\n' % (func_name, line_no))
        outfile.flush()
    
    def trace(frame, event, arg):
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
        for interesting in FUNC_TRACE_PATTERNS:
            if interesting in caller_file_name:
                want_to_know = True
        
        if want_to_know:
            # Must NEVER run when stdout.write() or flush() will trigger this trace
            outfile.write('Call to %s (%s @ %s) from %s (%s @ %s)\n'
                          % (func_name, file_name, line_no,
                             caller_name, caller_file_name, caller_line_no))
            outfile.flush()
            
            for interesting in LINE_TRACE_PATTERNS:
                if interesting in file_name:
                    return trace_lines
    
    sys.stdout.write('Enabling SUPER DEBUG - logging to OCELLARISSUPERDEBUG\n')
    sys.settrace(trace)
