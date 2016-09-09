"""
Plot timestep reports from an Ocellaris restart file 
"""
import h5py
import numpy
from matplotlib import pyplot
from matplotlib.widgets import Slider


def read_reports(h5_file_name, derived=True):
    hdf = h5py.File(h5_file_name, 'r')
    
    reps = {}
    for rep_name in hdf['/reports']:
        reps[rep_name] = numpy.array(hdf['/reports'][rep_name])
    
    if derived:
        if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
            reps['Et'] = reps['Ek'] + reps['Ep']  
    
    return reps


def plot_reports(file_name):
    reps = read_reports(file_name)
    report_names = sorted(reps.keys())
    N = len(report_names)
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    pyplot.subplots_adjust(left=0.15, bottom=0.25)
    
    line, = pyplot.plot([0], [0])
    
    ax_slider = fig.add_axes([0.15, 0.1, 0.65, 0.03])
    slider = Slider(ax_slider, 'Report', 0.5, N+0.499999, valinit=N/2)
    
    def update(val):
        i = int(round(val)-1)
        rep_name = report_names[i]
        x = reps['timesteps']
        y = reps[rep_name]
        
        line.set_data(x[-len(y):], y)
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
    h5_file_name = sys.argv[1]
    plot_reports(h5_file_name)
    pyplot.show()
