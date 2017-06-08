"""
Inspect timestep reports from one or more Ocellaris restart files 
"""
import os
import h5py
import numpy
import wx
import matplotlib
matplotlib.use('WxAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas, NavigationToolbar2WxAgg as NavigationToolbar


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
        if key == 'time':
            key = 'timesteps'
        reps[key] = arr
    
    if derived:
        if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
            N = min(reps['Ek'].size, reps['Ep'].size)
            reps['Et'] = reps['Ek'][:N] + reps['Ep'][:N]
    
    return reps


class OcellarisInspector(wx.Frame):
    def __init__(self, lables, report_names, reports):
        super(OcellarisInspector, self).__init__(None, title='Ocellaris Report Inspector')
        self.SetSize(800, 600)
        
        self.lables = lables
        self.report_names = report_names
        self.reports = reports
        
        self.layout_widgets()
        self.report_selected()

    def layout_widgets(self):
        p = wx.Panel(self)
        v = wx.BoxSizer(wx.VERTICAL)
        p.SetSizer(v)
        
        # Figure and figure controls
        self.fig = Figure((5.0, 4.0), dpi=100)
        self.canvas = FigureCanvas(p, wx.ID_ANY, self.fig)
        self.axes = self.fig.add_subplot(111)
        toolbar = NavigationToolbar(self.canvas)
        v.Add(self.canvas, proportion=1, flag=wx.EXPAND)
        v.Add(toolbar, flag=wx.EXPAND)
        v.Fit(self)
        
        # Choose report to plot
        h1 = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h1, flag=wx.ALL|wx.EXPAND, border=4)
        h1.Add(wx.StaticText(p, label='Report:'), flag=wx.CENTER)
        h1.AddSpacer(5)
        self.report_selector = wx.Choice(p, choices=self.report_names)
        self.report_selector.Select(0)
        self.report_selector.Bind(wx.EVT_CHOICE, self.report_selected)
        h1.Add(self.report_selector, proportion=1)
        
        # Customize the plot text
        h2 = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h2, flag=wx.ALL|wx.EXPAND, border=4)
        
        # Plot title
        h2.Add(wx.StaticText(p, label='Title:'), flag=wx.CENTER)
        h2.AddSpacer(5)
        self.title = wx.TextCtrl(p)
        self.title.Bind(wx.EVT_TEXT, self.update_plot)
        h2.Add(self.title, proportion=1)
        h2.AddSpacer(15)
        
        # Plot xlabel
        h2.Add(wx.StaticText(p, label='X:'), flag=wx.CENTER)
        h2.AddSpacer(5)
        self.xlabel = wx.TextCtrl(p)
        self.xlabel.Bind(wx.EVT_TEXT, self.update_plot)
        h2.Add(self.xlabel, proportion=1)
        h2.AddSpacer(15)
        
        # Plot ylabel
        h2.Add(wx.StaticText(p, label='Y:'), flag=wx.CENTER)
        h2.AddSpacer(5)
        self.ylabel = wx.TextCtrl(p)
        self.ylabel.Bind(wx.EVT_TEXT, self.update_plot)
        h2.Add(self.ylabel, proportion=1)
        h2.AddSpacer(5)
        
        
        v.Fit(p)
    
    def report_selected(self, evt=None):
        irep = self.report_selector.GetSelection()
        report_name = self.report_names[irep]
        
        self.title.ChangeValue('Ocellaris report %s' % report_name)
        self.xlabel.ChangeValue('t')
        self.ylabel.ChangeValue(report_name)
        
        self.update_plot()
        
    def update_plot(self, evt=None):
        irep = self.report_selector.GetSelection()
        report_name = self.report_names[irep]
        
        self.axes.clear()
        
        for i, label in enumerate(self.lables):
            x = self.reports[i]['timesteps']
            if report_name in self.reports[i]:
                y = self.reports[i][report_name]
            else:
                y = numpy.zeros_like(x)
                y[:] = numpy.NaN
            self.axes.plot(x, y, label=label)
        
        self.axes.relim()
        self.axes.autoscale_view()
        self.axes.set_title(self.title.GetValue())
        self.axes.set_xlabel(self.xlabel.GetValue())
        self.axes.set_ylabel(self.ylabel.GetValue())
        self.axes.legend(loc='best')
        self.fig.tight_layout()
        
        self.canvas.draw()


def show_inspector(file_names, lables):
    """
    Show wxPython window that allows chosing  which report to show
    """
    all_reps = []
    all_rep_names = set()
    for fn in file_names:
        reps = read_reports(fn)
        all_reps.append(reps)
        all_rep_names.update(reps.keys())
    report_names = sorted(all_rep_names)
    
    app = wx.App()
    frame = OcellarisInspector(lables, report_names, all_reps)
    frame.Show()
    app.MainLoop()


if __name__ == '__main__':
    import sys
    
    # Get report files to save
    h5_file_names = sys.argv[1:]
    
    # Get lables
    lables = []
    for i in range(len(h5_file_names)):
        fn = h5_file_names[i]
        if ':' in fn:
            fn, label = fn.split(':')
            h5_file_names[i] = fn
        else:
            bname = os.path.basename(fn)
            bname_split = bname.split('_endpoint_')
            label = bname_split[0]
        lables.append(label)
    
    # Make plots
    show_inspector(h5_file_names, lables)
