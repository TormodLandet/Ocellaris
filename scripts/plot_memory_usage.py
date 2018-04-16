"""
Plot the memory usage for an Ocellaris simulation based on log file data. You
must have specified output/show_memory_usage: yes in the input file to have
the MAX RSS memory information available 
"""
import sys
from collections import namedtuple
import numpy
from matplotlib import pyplot
from matplotlib.dates import date2num


MemoryLogStatement = namedtuple('MemoryLogStatement',
                                'timestep,iso_dt,max_rss,maj_faults,log_text')


def plot_memory_usage(log_statements):
    """
    Plot a list of MemoryLogStatement objects
    """
    max_rss = [mls.max_rss for mls in log_statements]
    maj_faults = [mls.maj_faults for mls in log_statements]
    datetimes = [numpy.datetime64(mls.iso_dt) for mls in log_statements]
    annotations = [mls.log_text for mls in log_statements] 
    
    max_rss = numpy.array(max_rss, dtype=float)
    maj_faults = numpy.array(maj_faults, dtype=float)
    datetimes = numpy.array(datetimes, dtype=numpy.datetime64).astype('O')
    
    ax = plot_annotated(datetimes, max_rss, annotations)
    ax.set_title('MAX RSS')
    ax.set_ylabel('MB')
    
    ax = plot_annotated(datetimes, maj_faults, annotations)
    ax.set_title('Major page faults')
    
    pyplot.show()


def plot_annotated(x, y, annotations):
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    line, = ax.plot(x, y, marker='x')
    
    # Information shown on hover
    info = ax.annotate("to be filled in", xy=(0,0), xytext=(-20,20), 
                       textcoords="offset points")
    info.set_visible(False)
    
    prev_i = -1
    def hover(event):
        nonlocal prev_i
        visible = info.get_visible()
        if event.inaxes == ax:
            contained, ind = line.contains(event)
            if contained:
                x, y = line.get_data()
                i = ind["ind"][0]
                if i == prev_i:
                    return
                prev_i = i 
                
                info.xy = (date2num(x[i]), y[i])
                info.set_text(annotations[i])
                info.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if visible:
                    info.set_visible(False)
                    fig.canvas.draw_idle()
    
    fig.canvas.mpl_connect("motion_notify_event", hover)
    return ax


def show_worst_offenders(log_statements):
    mls = log_statements[0]
    prev_max_rss = mls.max_rss
    print('-------------------------------------------------------------------')
    print('Start with MAX RSS %.1f MB at %s' % (mls.max_rss, mls.iso_dt))
    
    for mls in log_statements:
        if mls.max_rss > 10 + prev_max_rss:
            pst = (mls.max_rss - prev_max_rss) / prev_max_rss
            print('Increase from %.1f to %.1f (%.1f%%) at %s:\n\t%s'
                  % (prev_max_rss, mls.max_rss, pst * 100,
                     mls.iso_dt, mls.log_text))
        prev_max_rss = mls.max_rss
    
    print('End with MAX RSS %.1f MB at %s' % (mls.max_rss, mls.iso_dt))
    
    ############################################################################
    
    mls = log_statements[0]
    prev_maj_faults = mls.maj_faults
    print('-------------------------------------------------------------------')
    print('Start with maj. page faults %r at %s' % (mls.maj_faults, mls.iso_dt))
    
    for mls in log_statements:
        if mls.maj_faults > 5 + prev_maj_faults:
            pst = (mls.maj_faults - prev_maj_faults) / prev_maj_faults
            print('Increase from %r to %r (%.1f%%) at %s:\n\t%s'
                  % (prev_maj_faults, mls.maj_faults, pst * 100,
                     mls.iso_dt, mls.log_text))
        prev_maj_faults = mls.maj_faults
    
    print('End with maj. page faults %r at %s' % (mls.maj_faults, mls.iso_dt))
    print('-------------------------------------------------------------------')


def read_memory_usage(log_file_name, only_last_run=False):
    """
    Read an Ocellaris log file produced with 
    
        output:
            show_memory_usage: yes
    
    specified in the Ocellaris input file 
    """
    # Read the whole log as one large string
    log = open(log_file_name, 'rt').read()
    if only_last_run:
        log = log.split('Installed at:')[-1]
    
    # Split the log file string into records and store them
    raw_data = log.split('Current max RSS memory at timestep')
    log_statements = []
    text = raw_data[0].strip()
    for rd in raw_data[1:]:
        # Parse the memory log info
        words = rd[:100].split()
        timestep = int(words[0])
        iso_dt = words[2]
        max_rss = int(words[4])
        maj_faults = int(words[6])
        
        # Convert max RSS to MB. Assumes that ru_max_rss is in kB which is
        # correct on linux and BSD, but it is bytes on Mac OS X for some reason
        max_rss /= 1024
                
        # Store the log statement along with the preceding log text
        mls = MemoryLogStatement(timestep, iso_dt, max_rss, maj_faults, text)
        log_statements.append(mls)
        
        # The text of the next log statement
        i = rd.find('\n')
        text = '' if i == -1 else rd[i:].strip()
    
    return log_statements


if __name__ == '__main__':
    log_file_name = sys.argv[1]
    log_statements = read_memory_usage(log_file_name, only_last_run=True)
    show_worst_offenders(log_statements)
    plot_memory_usage(log_statements)
