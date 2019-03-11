# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
This plot panel is useful when a quasi static simulation contains multiple
quasi-static periods separated by transitions

An example is running a stair-stepped velocity and looking at the drag force
as it stabilises on each stair "thread" before moving on to the next velocity.

When comparing such simulations it is useful to syncronize the start of each
"stair" in case these do not happen at the same time instance for all
simulations
"""
import numpy
import matplotlib

matplotlib.use('WxAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigureCanvas,
    NavigationToolbar2WxAgg as NavigationToolbar,
)
import wx
from . import pub, TOPIC_METADATA, TOPIC_RELOAD
from .widget_shared import PlotLimSelectors, PlotCustomLine
from .panel_results import DEFAULT_REPORT, compute_new_lim
from ..algorithms import get_stairs


DEFAULT_CONTROL = 'Vel'


class OcellarisStairsPanel(wx.Panel):
    def __init__(self, parent, inspector_state):
        super().__init__(parent)
        self.istate = inspector_state
        self.need_update = True
        self.layout_widgets()
        self.reload_selector_data()

        # Select the default control parameter
        if DEFAULT_CONTROL in self.report_names:
            self.control_selector.Select(self.report_names.index(DEFAULT_CONTROL))
        else:
            self.control_selector.Select(0)

        # Select the default report to plot
        if DEFAULT_REPORT in self.report_names:
            self.report_selector.Select(self.report_names.index(DEFAULT_REPORT))
        else:
            self.report_selector.Select(0)

        self.control_selected()
        self.report_selected()

        self.Bind(wx.EVT_IDLE, self.on_idle)
        pub.subscribe(self.update_plot_soon, TOPIC_METADATA)
        pub.subscribe(self.reload_selector_data, TOPIC_RELOAD)

    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)

        # Figure and figure controls
        self.fig = Figure((5.0, 4.0), dpi=100)
        self.canvas = FigureCanvas(self, wx.ID_ANY, self.fig)
        self.axes = self.fig.add_subplot(111)
        toolbar = NavigationToolbar(self.canvas)
        self.plot_cursor_position_info = wx.StaticText(
            self, style=wx.ALIGN_RIGHT, size=(250, -1), label=''
        )
        self.canvas.mpl_connect('motion_notify_event', self.mouse_position_on_plot)
        v.Add(self.canvas, proportion=1, flag=wx.EXPAND)
        h = wx.BoxSizer(wx.HORIZONTAL)
        h.Add(toolbar, proportion=1)
        h.AddSpacer(10)
        h.Add(self.plot_cursor_position_info, flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(5)
        v.Add(h, flag=wx.EXPAND)

        # Choose report to plot
        gbs = wx.GridBagSizer(vgap=3, hgap=10)
        v.Add(gbs, flag=wx.ALL | wx.EXPAND, border=6)

        L = 0
        gbs.Add(
            wx.StaticText(self, label='Control report:'),
            flag=wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 0),
            span=(1, 1),
        )
        self.control_selector = wx.Choice(self)
        self.control_selector.Bind(wx.EVT_CHOICE, self.control_selected)
        gbs.Add(
            self.control_selector,
            flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 1),
            span=(1, 1),
        )

        gbs.Add(
            wx.StaticText(self, label='Control value:'),
            flag=wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 2),
            span=(1, 1),
        )
        self.control_value = wx.Choice(self)
        self.control_value.Bind(wx.EVT_CHOICE, self.control_value_selected)
        gbs.Add(
            self.control_value,
            flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 3),
            span=(1, 1),
        )

        L += 1
        gbs.Add(
            wx.StaticText(self, label='Plot report:'),
            flag=wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 0),
            span=(1, 1),
        )
        self.report_selector = wx.Choice(self)
        self.report_selector.Bind(wx.EVT_CHOICE, self.report_selected)
        gbs.Add(
            self.report_selector,
            flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 1),
            span=(1, 3),
        )

        # Plot title
        L += 1
        gbs.Add(
            wx.StaticText(self, label='Plot title:'),
            flag=wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 0),
            span=(1, 1),
        )
        self.title = wx.TextCtrl(self)
        self.title.Bind(wx.EVT_TEXT, self.update_plot_soon)
        gbs.Add(
            self.title,
            flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 1),
            span=(1, 3),
        )

        # Plot xlabel / log x axis
        L += 1
        gbs.Add(
            wx.StaticText(self, label='Label X:'),
            flag=wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 0),
            span=(1, 1),
        )
        self.xlabel = wx.TextCtrl(self)
        self.xlabel.Bind(wx.EVT_TEXT, self.update_plot_soon)
        gbs.Add(
            self.xlabel,
            flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 1),
            span=(1, 2),
        )
        self.xlog = wx.CheckBox(self, label='X as log axis')
        self.xlog.Bind(wx.EVT_CHECKBOX, self.update_plot_soon)
        gbs.Add(self.xlog, flag=wx.ALIGN_CENTER_VERTICAL, pos=(L, 3), span=(1, 1))

        # Plot ylabel
        L += 1
        gbs.Add(
            wx.StaticText(self, label='Label Y:'),
            flag=wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 0),
            span=(1, 1),
        )
        self.ylabel = wx.TextCtrl(self)
        self.ylabel.Bind(wx.EVT_TEXT, self.update_plot_soon)
        gbs.Add(
            self.ylabel,
            flag=wx.EXPAND | wx.ALIGN_CENTER_VERTICAL,
            pos=(L, 1),
            span=(1, 2),
        )
        self.ylog = wx.CheckBox(self, label='Y as log axis')
        self.ylog.Bind(wx.EVT_CHECKBOX, self.update_plot_soon)
        gbs.Add(self.ylog, flag=wx.ALIGN_CENTER_VERTICAL, pos=(L, 3), span=(1, 1))

        # Collapsible settings pane
        coll = wx.CollapsiblePane(self, label='Details')
        collp = coll.GetPane()
        coll.Bind(wx.EVT_COLLAPSIBLEPANE_CHANGED, lambda _: self.Layout())
        collv = wx.BoxSizer(wx.VERTICAL)
        collp.SetSizer(collv)
        v.Add(coll, flag=wx.EXPAND, proportion=0)

        self.plot_limits = PlotLimSelectors(collp, self.update_plot_soon)
        collv.Add(self.plot_limits, flag=wx.ALL | wx.EXPAND, border=4)

        self.custom_line = PlotCustomLine(collp, self.update_plot_soon)
        collv.Add(self.custom_line, flag=wx.ALL | wx.EXPAND, border=4)

        gbs.AddGrowableCol(1, proportion=1)
        v.Fit(self)

    def reload_selector_data(self, evt=None):
        # Store previous selection if found
        selected_ctrl = self.control_selector.GetSelection()
        if selected_ctrl != wx.NOT_FOUND:
            selected_ctrl = self.report_names[selected_ctrl]

        selected_rep = self.report_selector.GetSelection()
        if selected_rep != wx.NOT_FOUND:
            selected_rep = self.report_names[selected_rep]

        # Sort reports first by length then by name. This makes sure that
        # inner iteration reports end up last in the list
        all_rep_names = set()
        all_rep_lengths = {}
        for res in self.istate.results:
            rep_names = list(res.reports.keys())
            all_rep_names.update(rep_names)
            for rep_name in rep_names:
                all_rep_lengths[rep_name] = max(
                    all_rep_lengths.get(rep_name, 0), len(res.reports[rep_name])
                )

        def sort_key(rep_name):
            return (all_rep_lengths[rep_name], rep_name)

        self.report_names = sorted(all_rep_names, key=sort_key)

        # Update report selectors, restoring previous value if possible

        self.control_selector.Set(self.report_names)
        if selected_ctrl in self.report_names:
            self.control_selector.SetSelection(self.report_names.index(selected_ctrl))
        else:
            self.control_selector.SetSelection(0)
        self.control_selected()

        self.report_selector.Set(self.report_names)
        if selected_rep in self.report_names:
            self.report_selector.SetSelection(self.report_names.index(selected_rep))
        else:
            self.report_selector.SetSelection(0)
        self.report_selected()

        # In case the wx.Choice was empty it must now grow a bit vertically
        self.GetSizer().Layout()

    def mouse_position_on_plot(self, mpl_event):
        x, y = mpl_event.xdata, mpl_event.ydata
        if x is None or y is None:
            info = ''
        else:
            info = 'pos = (%.6g, %.6g)' % (x, y)

        self.plot_cursor_position_info.Label = info

    def control_selected(self, evt=None):
        irep = self.control_selector.GetSelection()
        if irep == wx.NOT_FOUND:
            self.control_values = []
            self.control_value.Set(self.control_values)
            self.control_value_selected()
            return
        report_name = self.report_names[irep]

        # Get previously selected value
        selected_val = self.control_value.GetSelection()
        if selected_val != wx.NOT_FOUND:
            selected_val = self.control_values[selected_val]

        # Update stair value selector
        values, limits = get_stairs(self.istate.results, report_name)
        self.control_values = [str(v) for v in values]
        self.stair_lims = limits
        self.control_value.Set(self.control_values)
        if selected_val in self.control_values:
            self.control_value.SetSelection(self.control_values.index(selected_val))
        else:
            self.control_value.SetSelection(0)

        self.control_value_selected()

    def control_value_selected(self, evt=None):
        ival = self.control_value.GetSelection()
        if ival == wx.NOT_FOUND:
            self.startstop_indices = None
        else:
            self.startstop_indices = self.stair_lims[ival]
        self.update_plot()

    def report_selected(self, evt=None):
        irep = self.report_selector.GetSelection()
        if irep == wx.NOT_FOUND:
            return
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
        self.axes.clear()

        if not self.startstop_indices:
            self.canvas.draw()
            self.need_update = False
            return

        irep = self.report_selector.GetSelection()
        if irep == wx.NOT_FOUND:
            self.canvas.draw()
            self.need_update = False
            return

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

        xs, ys = [], []
        for results, startstop in zip(
            self.istate.active_results, self.startstop_indices
        ):
            if report_name not in results.reports:
                plot([0], [None], label=results.label)
                continue
            istart_trans, istart, iend = startstop
            x = results.reports_x[report_name][istart_trans : iend + 1]
            x = numpy.array(x) - x[istart - istart_trans]
            y = results.reports[report_name][istart_trans : iend + 1]
            plot(x, y, label=results.label)
            xs.append(x)
            ys.append(y)

        self.axes.relim()
        self.axes.autoscale_view()

        if self.custom_line.active:
            xlim = self.axes.get_xlim()
            x = numpy.linspace(xlim[0], xlim[1], 1000)
            y, name = self.custom_line.get_function(x)
            plot(x, y, label=name, color='k', linestyle='--')
            xs.append(x)
            ys.append(y)

        # Change ylim if xlim is given
        xlim = self.plot_limits.get_xlim()
        ylim = self.plot_limits.get_ylim()
        if xlim != (None, None):
            xlocator = self.axes.xaxis.get_major_locator()
            ylim = compute_new_lim(xlim, ylim, xs, ys, xlocator)

        self.axes.set_xlim(xlim)
        self.axes.set_ylim(ylim)
        self.axes.set_title(self.title.GetValue())
        self.axes.set_xlabel(self.xlabel.GetValue())
        self.axes.set_ylabel(self.ylabel.GetValue())

        if len(self.istate.results) > 1:
            self.axes.legend(loc='best')
        self.fig.tight_layout()

        self.canvas.draw()
        self.need_update = False
