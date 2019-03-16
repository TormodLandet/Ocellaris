# Copyright (C) 2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

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


class OcellarisLineProbePanel(wx.Panel):
    def __init__(self, parent, inspector_state):
        super().__init__(parent)
        self.istate = inspector_state
        self.layout_widgets()
        self.reset_line_probes()

        self.Bind(wx.EVT_IDLE, self.on_idle)
        pub.subscribe(self.update_plot_soon, TOPIC_METADATA)
        pub.subscribe(self.reset_line_probes, TOPIC_RELOAD)

    def reset_line_probes(self, evt=None):
        if self.line_probe_selector.GetCount():
            self.line_probe_selector.Select(0)
            iprobe = self.line_probe_selector.GetSelection()
            selected_name = self.line_probe_names[iprobe]
        else:
            selected_name = ''

        all_probe_names = set()
        for res in self.istate.results:
            all_probe_names.update(res.line_probes.keys())
        self.line_probe_names = sorted(all_probe_names)

        idx = self.line_probe_selector.FindString(selected_name)
        if idx == wx.NOT_FOUND:
            idx = 0

        self.line_probe_selector.Set(self.line_probe_names)
        self.line_probe_selector.SetSelection(idx)

        self.active_line_probe = None
        self.need_update = True

        self.xlabel.ChangeValue('x-label')
        self.ylabel.ChangeValue('y-label')

        # In case the wx.Choice was empty it must now grow a bit vertically
        self.GetSizer().Layout()

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

        # Choose line probe to plot
        h = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h, flag=wx.ALL | wx.EXPAND, border=4)
        h.Add(wx.StaticText(self, label='Line probe:'), flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(5)
        self.line_probe_selector = wx.Choice(self)
        self.line_probe_selector.Bind(wx.EVT_CHOICE, self.on_select_line_probe)
        h.Add(self.line_probe_selector, proportion=1)
        h.AddSpacer(5)
        b = wx.Button(self, label='Load')
        b.Bind(wx.EVT_BUTTON, self.load_line_probes)
        h.Add(b)

        # Choose timestep
        h = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h, flag=wx.ALL | wx.EXPAND, border=4)
        h.Add(wx.StaticText(self, label='Time step:'), flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(5)
        self.time_selector = wx.Slider(self)

        def lam(evt):
            return setattr(self, 'need_update', True)

        self.time_selector.Bind(wx.EVT_SLIDER, lam)
        h.Add(self.time_selector, proportion=1)

        # Customize the plot text
        fgs = wx.FlexGridSizer(cols=2, vgap=3, hgap=10)
        fgs.AddGrowableCol(1, proportion=1)
        v.Add(fgs, flag=wx.ALL | wx.EXPAND, border=6)

        # Plot title
        fgs.Add(wx.StaticText(self, label='Plot title:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.title = wx.TextCtrl(self)
        self.title.Bind(wx.EVT_TEXT, self.update_plot_soon)
        fgs.Add(self.title, flag=wx.EXPAND)

        # Plot xlabel / log x axis
        fgs.Add(wx.StaticText(self, label='Label X:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.xlabel = wx.TextCtrl(self)
        self.xlabel.Bind(wx.EVT_TEXT, self.update_plot_soon)
        fgs.Add(self.xlabel, flag=wx.EXPAND)

        # Plot ylabel
        fgs.Add(wx.StaticText(self, label='Label Y:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ylabel = wx.TextCtrl(self)
        self.ylabel.Bind(wx.EVT_TEXT, self.update_plot_soon)
        fgs.Add(self.ylabel, flag=wx.EXPAND)

        # Plot ylabel
        fgs.Add(wx.StaticText(self, label='Aspect ratio:'), flag=wx.ALIGN_CENTER_VERTICAL)
        self.aspect = wx.Choice(self)
        self.aspect.Set(['equal', 'auto'])
        self.aspect.Select(1)
        self.aspect.Bind(wx.EVT_CHOICE, self.update_plot_soon)
        fgs.Add(self.aspect, flag=wx.EXPAND)

        v.Fit(self)

    def mouse_position_on_plot(self, mpl_event):
        x, y = mpl_event.xdata, mpl_event.ydata
        if x is None or y is None:
            info = ''
        else:
            info = 'pos = (%.6g, %.6g)' % (x, y)
        self.plot_cursor_position_info.Label = info

    def load_line_probes(self, evt=None):
        iprobe = self.line_probe_selector.GetSelection()
        line_probe_name = self.line_probe_names[iprobe]
        self.active_line_probe = line_probe_name

        max_num_timesteps = 0
        tmin, tmax = 1e100, -1e100

        # Populate the cache for these line probes
        for results in self.istate.active_results:
            if line_probe_name in results.line_probes:
                probe = results.line_probes[line_probe_name]

                with wx.BusyCursor():
                    try:
                        _description, timesteps, _xvec, _yvec, _zvec, _data = probe.get_line_probes(
                            cache=True
                        )
                    except BaseException:
                        continue

                    if len(timesteps) == 0:
                        continue

                    # Compute bounds
                    tmin = min(tmin, timesteps[0])
                    tmax = max(tmax, timesteps[-1])
                    max_num_timesteps = max(max_num_timesteps, len(timesteps))

        self.tmin = tmin
        self.tmax = tmax
        self.max_num_timesteps = max_num_timesteps
        self.time_selector.SetMax(self.max_num_timesteps)

        self.title.ChangeValue('Line probe %s' % line_probe_name)
        self.xlabel.ChangeValue('x-label')
        self.ylabel.ChangeValue('y-label')

        self.update_plot()

    def on_select_line_probe(self, evt=None):
        self.active_line_probe = None
        self.need_update = True

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
        if self.active_line_probe is None or self.max_num_timesteps == 0:
            self.canvas.draw()
            return
        line_probe_name = self.active_line_probe

        # Get time
        f = self.time_selector.GetValue() * 1.0 / self.max_num_timesteps
        t = self.tmin + (self.tmax - self.tmin) * f

        tsel = set()
        abcissa_label = ''
        plotted_target = False
        for ires, results in enumerate(self.istate.active_results):
            if line_probe_name in results.line_probes:
                probe = results.line_probes[line_probe_name]
                with wx.BusyCursor():
                    try:
                        _description, timesteps, xvec, yvec, zvec, data = probe.get_line_probes(
                            cache=True
                        )
                    except Exception:
                        continue

                i = numpy.argmin(abs(timesteps - t))
                dt = timesteps[1] - timesteps[0]
                if abs(timesteps[i] - t) > 1.5 * dt:
                    continue
                tsel.add(timesteps[i])

                # For plotting, figure out which axis is the abcissa
                if xvec[0] != xvec[-1]:
                    abcissa = xvec
                    abcissa_label = 'x-pos'
                elif yvec[0] != yvec[-1]:
                    abcissa = yvec
                    abcissa_label = 'y-pos'
                else:
                    abcissa = zvec
                    abcissa_label = 'z-pos'

                self.axes.plot(abcissa, data[i], color='C%d' % ires, label=results.label)

                if probe.target_ordinate and not plotted_target:
                    assert len(probe.target_ordinate) == len(probe.target_abcissa)
                    self.axes.plot(
                        probe.target_abcissa, probe.target_ordinate, 'k:', label=probe.target_name
                    )
                    plotted_target = True

        self.axes.set_aspect(self.aspect.GetStringSelection())

        if len(tsel) == 1:
            t = tsel.pop()

        if self.xlabel.Value == 'x-label' and abcissa_label:
            self.xlabel.ChangeValue(abcissa_label)
        if self.ylabel.Value == 'y-label':
            self.ylabel.ChangeValue(line_probe_name)

        self.axes.set_title(self.title.GetValue())
        self.axes.set_xlabel(self.xlabel.GetValue())
        self.axes.set_ylabel(self.ylabel.GetValue())
        if len(self.istate.results) > 1 or plotted_target:
            self.axes.legend(loc='best')
        self.fig.tight_layout()

        self.canvas.draw()
        self.need_update = False
