"""
Plot convergence summary files 
"""
import json
import codecs
import collections
import math
from matplotlib import pyplot


def read_summary_file(summary_file_name,):
    summaries = collections.OrderedDict()
    triangles = []
    global_meta = {}
    
    inp = codecs.open(summary_file_name, 'r', 'utf-8')
    line = inp.readline()
    while line:
        if line.startswith('#' * 40):
            # This is the start of a header, the header is the next line
            name = inp.readline().strip()
            inp.readline() # read the line after the header
            summaries[name] = discr, errL2, errH1, meta, style = [], [], [], {}, {}
            meta['y-label'] = 'L2 error'
            style['label'] = name
            style['marker'] = '.'
            L2 = True
            
        elif line.startswith('{'):
            s = json.loads(line)
            style.update(s)
        
        elif line.startswith('GLOBAL'):
            gm = json.loads(line[7:])
            global_meta.update(gm)
        
        elif line.startswith('HIDE'):
            del summaries[name]
        
        elif line.startswith('TRIANGLE'):
            wds = line.split()
            trislope = float(wds[1])
            x0, y0, x1 = [float(w) for w in wds[2:8]]
            y2 = y0*math.exp(trislope*math.log(x0/x1))
            xt = [x0, x1, x0, x0]
            yt = [y0, y0, y2, y0]
            
            ts2 = math.log(y0/y2)/math.log(x0/x1) 
            print 'Triangle:', trislope, x0, y0, x1, ts2
            triangles.append(('Order %g' % trislope, xt, yt))
        
        elif line.startswith(' Discr.'):
            # We are at the table of results
            inp.readline()
            discr_type = inp.readline().split()[0]
            meta['x-label'] = discr_type
            inp.readline()
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
    
    return summaries, triangles, global_meta


def plot_convergence_summaries(file_name):
    summaries, triangles, global_meta = read_summary_file(file_name)
    case_names = summaries.keys()
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Comparison of convergence tests')
    
    # Plot the convergence plot series in loglog
    for name in case_names:
        discr, errL2, errH1, meta, style = summaries[name]
        
        print name
        print discr, errL2, errH1
        print style
        print
        
        ax.loglog(discr, errL2, **style)
        
    # Plot any convergence order triangles
    for triname, xt, yt in triangles:
        ax.loglog(xt, yt, '-k', label=None, lw=1)
        ax.text(xt[0], yt[0], triname, va='bottom')
    
    ax.set_xlabel(meta['x-label'])
    ax.set_ylabel(meta['y-label'])
    ax.legend()
    
    for attr, val in global_meta.items():
        getattr(ax, 'set_' + attr)(val)
    
    fig.tight_layout()


if __name__ == '__main__':
    import sys
    summary_file_name = sys.argv[1]
    plot_convergence_summaries(summary_file_name)
    pyplot.show()
