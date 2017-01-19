"""
Plot convergence summary files 
"""
import json
import codecs
import collections
from matplotlib import pyplot


def read_summary_file(summary_file_name,):
    summaries = collections.OrderedDict()
    
    inp = codecs.open(summary_file_name, 'r', 'utf-8')
    line = inp.readline()
    while line:
        if line.startswith('#'*40):
            # This is the start of a header, the header is the next line
            name = inp.readline().strip()
            inp.readline() # read the line after the header
            summaries[name] = discr, errL2, errH1, style = [], [], [], {}
            style['label'] = name
            L2 = True
            
        elif line.startswith('{'):
            s = json.loads(line)
            style.update(s)
        
        elif line.startswith(' Discr.'):
            # We are at the table of results
            inp.readline(); inp.readline(); inp.readline()
            while True:
                wds = inp.readline().split()
                try:
                    # Read dicr and u0 + u1 columns
                    d, e = float(wds[0]), float(wds[1]) + float(wds[2])
                    if L2:
                        discr.append(d)
                        errL2.append(e)
                    else:
                        errH1.append(e)
                except ValueError:
                    break
            L2 = False
        
        line = inp.readline()
    
    return summaries


def plot_convergence_summaries(file_name):
    summaries = read_summary_file(file_name)
    case_names = summaries.keys()
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Comparison of convergence tests')
    
    for name in case_names:
        discr, errL2, errH1, style = summaries[name]
        
        print name
        print discr, errL2, errH1
        print style
        print
        
        ax.loglog(discr, errL2, **style)
    
    ax.legend()
    fig.tight_layout()


if __name__ == '__main__':
    import sys
    summary_file_name = sys.argv[1]
    plot_convergence_summaries(summary_file_name)
    pyplot.show()
