import sys

def report_error(header, description, stop=True):
    print 'ERROR === '*8
    print
    print header
    print
    print description
    print
    print 'ERROR === '*8 
    if stop:
        print 'Stopping due to error'
        sys.exit(-1)
    