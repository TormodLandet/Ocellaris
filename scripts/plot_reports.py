"""
Plot timestep reports from one or more Ocellaris restart files 
"""
import os
import h5py
import numpy
from matplotlib import pyplot
from matplotlib.widgets import Slider


def read_reports(file_name, derived=True):
    if file_name.endswith('h5'):
        return read_reports_h5(file_name, derived)
    else:
        return read_reports_log(file_name, derived)


def read_reports_h5(h5_file_name, derived=True):
    hdf = h5py.File(h5_file_name, 'r')
    
    reps = {}
    for rep_name in hdf['/reports']:
        reps[rep_name] = numpy.array(hdf['/reports'][rep_name])
    
    if derived:
        if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
            reps['Et'] = reps['Ek'] + reps['Ep']  
    
    return reps


def read_reports_log(log_file_name, derived=True):
    data = {}
    for line in open(log_file_name, 'rt'):
        if line.startswith('Reports for timestep'):
            parts = line[12:].split(',')
            for pair in parts:
                try:
                    key, value = pair.split('=')
                    key = key.strip()
                    value = float(value)
                    data.setdefault(key, []).append(value)
                except:
                    break
    
    reps = {}
    for key, values in data.items():
        arr = numpy.array(values)
        if key == 'timestep':
            key = 'timesteps'
        reps[key] = arr
    
    if derived:
        if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
            reps['Et'] = reps['Ek'] + reps['Ep']  
    
    return reps


def plot_reports(file_names):
    all_reps = []
    all_rep_names = set()
    for fn in file_names:
        reps = read_reports(fn)
        all_reps.append(reps)
        all_rep_names.update(reps.keys())
    Nfiles = len(file_names)
    
    report_names = sorted(all_rep_names)
    N = len(report_names)
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    pyplot.subplots_adjust(left=0.15, bottom=0.25)
    
    lines = []
    for fn in file_names:
        bname = os.path.basename(fn)
        bname_split = bname.split('_endpoint_')
        label = bname_split[0]
        line, = pyplot.plot([0], [0], label=label)
        lines.append(line)
    if len(file_names) > 1:
        ax.legend()
    
    ax_slider = fig.add_axes([0.15, 0.1, 0.65, 0.03])
    slider = Slider(ax_slider, 'Report', 0.5, N+0.499999, valinit=N/2)
    
    def update(val):
        i = int(round(val)-1)
        rep_name = report_names[i]
        
        for i in range(Nfiles):
            x = all_reps[i]['timesteps']
            if rep_name in all_reps[i]:
                y = all_reps[i][rep_name]
            else:
                y = numpy.zeros_like(x)
                y[:] = numpy.NaN
            lines[i].set_data(x[-len(y):], y)
        
        ax.relim()
        ax.autoscale_view()
        ax.set_title('Ocellaris report %s' % rep_name)
        ax.set_ylabel(rep_name)
        slider.valtext.set_text(rep_name)
        
        fig.canvas.draw()
    
    slider.on_changed(update)
    update(slider.val)
    
    ax.set_xlabel('time')
    
    # avoid garbage collection of the slider making it non responsive
    return slider


if __name__ == '__main__':
    import sys
    h5_file_names = sys.argv[1:]
    plot_reports(h5_file_names)
    pyplot.show()
