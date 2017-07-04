"""
Inspect timestep reports from one or more Ocellaris restart files 
"""
import numpy
import wx
from wx.lib.pubsub import pub
import matplotlib
matplotlib.use('WxAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas, NavigationToolbar2WxAgg as NavigationToolbar
from ocellaris_post import Results
from .icons import OCELLARIS_ICON


TOPIC_METADATA = 'metadata_updated'


class OcellarisInspector(wx.Frame):
    def __init__(self, results):
        super(OcellarisInspector, self).__init__(None, title='Ocellaris Report Inspector')
        self.results = results
        
        self.layout_widgets()
        self.SetSize(800, 800)
        self.SetIcon(OCELLARIS_ICON.GetIcon())
        
    def layout_widgets(self):
        p = wx.Panel(self)
        nb = wx.Notebook(p, style=wx.NB_BOTTOM)
        
        self.metadata_panel = OcellarisMetadataPanel(nb, self.results)
        nb.AddPage(self.metadata_panel, 'Metadata')
        self.metadata_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        self.reports_panel = OcellarisReportsPanel(nb, self.results)
        nb.AddPage(self.reports_panel, 'Timestep reports')
        self.reports_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        nb.SetSelection(1)
        s = wx.BoxSizer()
        s.Add(nb, 1, wx.EXPAND)
        p.SetSizer(s)


class OcellarisReportsPanel(wx.Panel):
    def __init__(self, parent, results):
        super(OcellarisReportsPanel, self).__init__(parent)
        self.results = results
        self.reports = []
        all_rep_names = set()
        for res in results:
            self.reports.append(res.reports)
            all_rep_names.update(res.reports.keys())
        self.report_names = sorted(all_rep_names)
        
        self.layout_widgets()
        self.report_selected()
        
        self.Bind(wx.EVT_IDLE, self.on_idle)
        pub.subscribe(self.update_plot_soon, TOPIC_METADATA)
    
    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)
        
        # Figure and figure controls
        self.fig = Figure((5.0, 4.0), dpi=100)
        self.canvas = FigureCanvas(self, wx.ID_ANY, self.fig)
        self.axes = self.fig.add_subplot(111)
        toolbar = NavigationToolbar(self.canvas)
        self.plot_cursor_position_info = wx.StaticText(self, style=wx.ALIGN_RIGHT, size=(250, -1), label='')
        self.canvas.mpl_connect('motion_notify_event', self.mouse_position_on_plot)
        v.Add(self.canvas, proportion=1, flag=wx.EXPAND)
        h = wx.BoxSizer(wx.HORIZONTAL)
        h.Add(toolbar, proportion=1)
        h.AddSpacer(10)
        h.Add(self.plot_cursor_position_info, flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(5)
        v.Add(h, flag=wx.EXPAND)
        
        # Choose report to plot
        h1 = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h1, flag=wx.ALL|wx.EXPAND, border=4)
        h1.Add(wx.StaticText(self, label='Report:'), flag=wx.ALIGN_CENTER_VERTICAL)
        h1.AddSpacer(5)
        self.report_selector = wx.Choice(self, choices=self.report_names)
        self.report_selector.Select(0)
        self.report_selector.Bind(wx.EVT_CHOICE, self.report_selected)
        h1.Add(self.report_selector, proportion=1)
        
        # Customize the plot text
        Nrows = len(self.results) + 3
        #Nrows2 = Nrows // 2 + 1 if Nrows % 2 else Nrows // 2
        fgs = wx.FlexGridSizer(rows=Nrows, cols=3, vgap=3, hgap=10)
        fgs.AddGrowableCol(1, proportion=1)
        #fgs.AddGrowableCol(4, proportion=1)
        v.Add(fgs, flag=wx.ALL|wx.EXPAND, border=6)
        
        # Plot title
        fgs.Add(wx.StaticText(self, label='Plot title:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.title = wx.TextCtrl(self)
        self.title.Bind(wx.EVT_TEXT, self.update_plot_soon)
        fgs.Add(self.title, flag=wx.EXPAND)
        fgs.AddSpacer(0)
        
        # Plot xlabel / log x axis
        fgs.Add(wx.StaticText(self, label='Label X:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.xlabel = wx.TextCtrl(self)
        self.xlabel.Bind(wx.EVT_TEXT, self.update_plot_soon)
        fgs.Add(self.xlabel, flag=wx.EXPAND)
        self.xlog = wx.CheckBox(self, label='X as log axis')
        self.xlog.Bind(wx.EVT_CHECKBOX, self.update_plot_soon)
        fgs.Add(self.xlog)
        
        # Plot ylabel
        fgs.Add(wx.StaticText(self, label='Label Y:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ylabel = wx.TextCtrl(self)
        self.ylabel.Bind(wx.EVT_TEXT, self.update_plot_soon)
        fgs.Add(self.ylabel, flag=wx.EXPAND)
        self.ylog = wx.CheckBox(self, label='Y as log axis')
        self.ylog.Bind(wx.EVT_CHECKBOX, self.update_plot_soon)
        fgs.Add(self.ylog)
        
        v.Fit(self)
        
    def mouse_position_on_plot(self, mpl_event):
        x, y = mpl_event.xdata, mpl_event.ydata
        if x is None or y is None:
            info = ''
        else:
            info = 'pos = (%.6g, %.6g)' % (x, y)
        self.plot_cursor_position_info.Label = info

    def report_selected(self, evt=None):
        irep = self.report_selector.GetSelection()
        report_name = self.report_names[irep]
        
        self.title.ChangeValue('Ocellaris report %s' % report_name)
        self.xlabel.ChangeValue('t')
        self.ylabel.ChangeValue(report_name)
        
        self.update_plot()
    
    def update_plot_soon(self, evt=None):
        """
        Update the plot the next time the event queue is empty
        """
        self.need_update = True
    
    def on_idle(self, evt=None):
        if self.need_update:
            self.update_plot()
    
        
    def update_plot(self, evt=None):
        """
        Update the plot at once
        """
        irep = self.report_selector.GetSelection()
        report_name = self.report_names[irep]
        
        # How to plot
        xlog = self.xlog.GetValue()
        ylog = self.ylog.GetValue()
        if xlog and ylog:
            plot = self.axes.loglog
        elif xlog:
            plot = self.axes.semilogx
        elif ylog:
            plot = self.axes.semilogy
        else:
            plot = self.axes.plot
        
        self.axes.clear()
        
        for i, results in enumerate(self.results):
            x = self.reports[i]['timesteps']
            if report_name in self.reports[i]:
                y = self.reports[i][report_name]
            else:
                y = numpy.zeros_like(x)
                y[:] = numpy.NaN
            
            plot(x, y, label=results.label)
        
        self.axes.relim()
        self.axes.autoscale_view()
        self.axes.set_title(self.title.GetValue())
        self.axes.set_xlabel(self.xlabel.GetValue())
        self.axes.set_ylabel(self.ylabel.GetValue())
        self.axes.legend(loc='best')
        self.fig.tight_layout()
        
        self.canvas.draw()
        self.need_update = False


class OcellarisMetadataPanel(wx.Panel):
    def __init__(self, parent, results):
        super(OcellarisMetadataPanel, self).__init__(parent)
        self.results = results
        self.layout_widgets()
    
    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)
        
        # Metadata FlexGridSizer
        Nrows = len(self.results)
        fgs = wx.FlexGridSizer(rows=Nrows, cols=3, vgap=3, hgap=10)
        fgs.AddGrowableCol(1, proportion=1)
        v.Add(fgs, flag=wx.ALL|wx.EXPAND, border=6)
        
        # Customize the lables
        self.label_controls = []
        for il, results in enumerate(self.results):
            fgs.Add(wx.StaticText(self, label='File %d label:' % il), flag=wx.ALIGN_CENTER_VERTICAL)
            label_ctrl = wx.TextCtrl(self, value=results.label)
            label_ctrl.Bind(wx.EVT_TEXT, self.update_lables)
            fgs.Add(label_ctrl, flag=wx.EXPAND)
            st = wx.StaticText(self, label='(%s)' % results.label)
            st.SetToolTip(results.file_name)
            fgs.Add(st, flag=wx.ALIGN_CENTER_VERTICAL)
            self.label_controls.append(label_ctrl)
        
        v.Fit(self)
    
    def update_lables(self, event):
        for label_control, results in zip(self.label_controls, self.results):
            results.label = label_control.GetValue()
        pub.sendMessage(TOPIC_METADATA)


def show_inspector(file_names, lables):
    """
    Show wxPython window that allows chosing  which report to show
    """
    results = []
    for file_name, label in zip(file_names, lables):
        res = Results(file_name)
        res.label = label
        results.append(res)
    
    app = wx.App()
    frame = OcellarisInspector(results)
    frame.Show()
    app.MainLoop()
