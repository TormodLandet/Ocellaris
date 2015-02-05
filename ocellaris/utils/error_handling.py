import sys

def report_error(header, description, stop=False):
    print 'ERROR === '*8
    print
    print header
    print
    print description
    print
    print 'ERROR === '*8 
    if stop:
        sys.exit(-1)
    