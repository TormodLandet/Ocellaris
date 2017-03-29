"""
Plot timestep reports from one or more Ocellaris restart files 
"""
import os
import StringIO
import urllib, base64
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
            N = min(reps['Ek'].size, reps['Ep'].size)
            reps['Et'] = reps['Ek'][:N] + reps['Ep'][:N]
    
    return reps


def plot_reports(file_names, save=False):
    """
    Show matplotlib window with a slider that allows chosing 
    which report to show. If save==True then save reports to
    png+html instead
    """
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
    
    
    
    lines = []
    for fn in file_names:
        bname = os.path.basename(fn)
        bname_split = bname.split('_endpoint_')
        label = bname_split[0]
        line, = pyplot.plot([0], [0], label=label)
        lines.append(line)
    if len(file_names) > 1:
        ax.legend()
    
    def plot_rep(rep_name):
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
        
        fig.canvas.draw()
    
    def update(val):
        i = int(round(val)-1)
        rep_name = report_names[i]
        plot_rep(rep_name)
        slider.valtext.set_text(rep_name)
    
    if save:
        # Save report instead of showing plots on screen
        html_file_name = 'reports.html'
        with open(html_file_name, 'wt') as html:
            html.write('<!DOCTYPE html>\n<html lang="en">\n<head>\n')
            html.write('  <meta charset="utf-8">\n')
            html.write('  <title>Ocellaris reports</title>\n')
            html.write('  <style>\n')
            html.write('    body * { margin-left: 5%; }\n')
            html.write('    img { margin-left: 5%; }\n')
            html.write('  </style>\n')
            html.write('</head>\n<body>\n')
            html.write('\n\n<h1>Ocellaris reports</h1>\n\n')
  
            fig.size = (8, 6)
            for rep_name in report_names:
                html.write('\n<h2>%s</h2>\n' % rep_name)
                plot_rep(rep_name)
                fig.tight_layout()
                
                # Get png data as base64 encoded <img> element
                imgdata = StringIO.StringIO()
                fig.savefig(imgdata, format='png')
                imgdata.seek(0)  # rewind the data
                png = base64.b64encode(imgdata.buf)
                html.write('<img alt="Ocellaris report %s" ' % rep_name +
                           'src="data:image/png;base64,%s">\n' % urllib.quote(png))
            html.write('\n\n</body>\n</html>')
            print 'Wrote report file', html_file_name
        return
    
    pyplot.subplots_adjust(left=0.15, bottom=0.25)
    ax_slider = fig.add_axes([0.15, 0.1, 0.65, 0.03])
    slider = Slider(ax_slider, 'Report', 0.5, N+0.499999, valinit=N/2)
    
    slider.on_changed(update)
    update(slider.val)
    
    ax.set_xlabel('time')
    
    # avoid garbage collection of the slider making it non responsive
    return slider


if __name__ == '__main__':
    import sys
    
    # Get report files to save
    h5_file_names = sys.argv[1:]
    
    save = False
    if '--save' in h5_file_names:
        h5_file_names.remove('--save')
        save = True
    
    plot_reports(h5_file_names, save)
    
    if not save:
        pyplot.show()
