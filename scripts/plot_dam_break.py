# encoding: utf-8
from __future__ import division
import os
import numpy
from matplotlib import pyplot


def read_iso_surface_file(file_name):
    with open(file_name, 'rt') as f:
        # Read header
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()

        field_name = line1[31:-7]
        value = float(line2[9:])
        dim = int(line3[8:])
        
        assert dim == 2
        times = []
        lines = []
        
        tline = f.readline()
        while tline:
            wt = tline.split()
            time = float(wt[1])
            nlines = int(wt[3])
            
            tlines = []
            for _ in range(nlines):
                xvals = [float(v) for v in f.readline().split()]
                yvals = [float(v) for v in f.readline().split()]
                zvals = [float(v) for v in f.readline().split()]
                tlines.append((xvals, yvals, zvals))
            
            times.append(time)
            lines.append(tlines)
            tline = f.readline()
        
        return field_name, value, dim, times, lines
    

def get_simulation_name(file_name, names_allready_taken):
    # Find a name for this simulation
    basename = os.path.basename(file_name)
    if '_free_surface' in basename:
        basename = basename.split('_free_surface')[0]
    i = 1
    name = basename
    while name in names_allready_taken:
        i += 1
        name = basename + str(i)
    names_allready_taken.add(name)
    return name


def plot_iso_surface_file(file_names, n=2**0.5, a=0.05715, g=9.81):
    """
    Plot free surface elevations in the format of Martin and Moyce (1952)
    The definition of the parameters n and a are from the same article
    and relate to the width and height of the fluid column
    
    :param list[str] file_names: names of Ocellaris iso surface ouput files 
    :param float n: height_column == n**2 * a
    :param float a: width_column
    """
    # Data from Martin and Moyce (1952), Table 2 and 6
    mmTvec = [0.41, 0.84, 1.19, 1.43, 1.63, 1.83, 1.98, 2.20, 2.32, 2.51, 2.65,
              2.83, 2.98, 3.11, 3.33]
    mmZvec = [1.11, 1.22, 1.44, 1.67, 1.89, 2.11, 2.33, 2.56, 2.78, 3.00, 3.22,
              3.44, 3.67, 3.89, 4.11]
    mmYvec = [0.56, 0.77, 0.93, 1.08, 1.28, 1.46, 1.66, 1.84, 2.00, 2.21, 2.45,
              2.70, 3.06, 3.44, 4.20, 5.25, 7.40]
    mmHvec = [0.94, 0.89, 0.83, 0.78, 0.72, 0.67, 0.61, 0.56, 0.50, 0.44, 0.39,
              0.33, 0.28, 0.22, 0.17, 0.11, 0.06]
    
    plots = [('Horizontal maximum', [('MM', mmTvec, mmZvec)]),
             ('Vertical maximum',   [('MM', mmYvec, mmHvec)])]
    
    # Read files
    names = set()
    for file_name in file_names:
        field_name, value, dim, times, lines = read_iso_surface_file(file_name)
        label = get_simulation_name(file_name, names) 
        print label
        print field_name, value, dim
        
        # Y = τ
        Tvec, Yvec, Zvec, Hvec = [], [], [], []
        for i, tlines in enumerate(lines):
            txmax = tymax = -numpy.inf
            for xvals, yvals, _zvals in tlines:
                if len(xvals): txmax = max(txmax, numpy.max(xvals))
                if len(yvals): tymax = max(tymax, numpy.max(yvals))
            Tvec.append(times[i]*(g/a)**0.5*n)
            Yvec.append(times[i]*(g/a)**0.5)
            Zvec.append(txmax/a)
            Hvec.append(tymax/(a*n**2))
        print 'tmax, Tmax, Ymax', times[-1], Tvec[-1], Yvec[-1]
        
        plots[0][1].append((label, Tvec, Zvec))
        plots[1][1].append((label, Yvec, Hvec))
    
    # Plot surface elevations with time
    for name, lines in plots:
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_title(name)
        
        for label, tvec, vals in lines:
            kwargs = dict(marker='o', ls='') if label == 'MM' else {}
            ax.plot(tvec, vals, label=label, **kwargs)
        if len(lines) > 1:
            ax.legend()


if __name__ == '__main__':
    import sys
    iso_surface_file_names = sys.argv[1:]
    plot_iso_surface_file(iso_surface_file_names)
    pyplot.show()
