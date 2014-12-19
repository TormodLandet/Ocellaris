import sys

def report_error(header, description, stop=False):
    print '='*80
    print 'ERROR:', header
    print
    print description
    print
    print '='*80 
    if stop:
        sys.exit(-1)
    