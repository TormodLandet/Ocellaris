import matplotlib
matplotlib.use('WxAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas, NavigationToolbar2WxAgg as NavigationToolbar
import wx
from . import pub, TOPIC_METADATA


DEFAULT_REPORT = 'Cof_max'


class OcellarisReportsPanel(wx.Panel):
    def __init__(self, parent, results):
        super(OcellarisReportsPanel, self).__init__(parent)
        self.results = results
        
        # Sort reports first by length then by name. This makes sure that 
        # inner iteration reports end up last in the list
        all_rep_names = set()
        all_rep_lengths = {}
        for res in results:
            rep_names = res.reports.keys()
            all_rep_names.update(rep_names)
            for rep_name in rep_names:
                all_rep_lengths[rep_name] = max(all_rep_lengths.get(rep_name, 0),
                                                len(res.reports[rep_name]))
        sort_key = lambda rep_name: (all_rep_lengths[rep_name], rep_name)
        self.report_names = sorted(all_rep_names, key=sort_key)
        self.need_update = True
        
        self.layout_widgets()
        
        # Select the default report
        if DEFAULT_REPORT in self.report_names:
            self.report_selector.Select(self.report_names.index(DEFAULT_REPORT))
        else:
            self.report_selector.Select(0)
        
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
        self.report_selector.Bind(wx.EVT_CHOICE, self.report_selected)
        h1.Add(self.report_selector, proportion=1)
        
        # Customize the plot text
        fgs = wx.FlexGridSizer(rows=3, cols=3, vgap=3, hgap=10)
        fgs.AddGrowableCol(1, proportion=1)
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
        
        for results in self.results:
            if report_name not in results.reports:
                plot([0], [None], label=results.label)
                continue
            x = results.reports_x[report_name]
            y = results.reports[report_name]
            plot(x, y, label=results.label)
                
        self.axes.relim()
        self.axes.autoscale_view()
        self.axes.set_title(self.title.GetValue())
        self.axes.set_xlabel(self.xlabel.GetValue())
        self.axes.set_ylabel(self.ylabel.GetValue())
        if len(self.results) > 1:
            self.axes.legend(loc='best')
        self.fig.tight_layout()
        
        self.canvas.draw()
        self.need_update = False
