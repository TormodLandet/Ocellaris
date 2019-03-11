# Copyright (C) 2017-2019 Tormod Landet
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


class OcellarisSurfacesPanel(wx.Panel):
    def __init__(self, parent, inspector_state):
        super().__init__(parent)
        self.istate = inspector_state
        self.layout_widgets()
        self.reset_surfaces()

        self.Bind(wx.EVT_IDLE, self.on_idle)
        pub.subscribe(self.update_plot_soon, TOPIC_METADATA)
        pub.subscribe(self.reset_surfaces, TOPIC_RELOAD)

    def reset_surfaces(self, evt=None):
        if self.surface_selector.GetCount():
            self.surface_selector.Select(0)
            isurf = self.surface_selector.GetSelection()
            selected_name = self.surface_names[isurf]
        else:
            selected_name = ''

        all_surface_names = set()
        for res in self.istate.results:
            all_surface_names.update(res.surfaces.keys())
        self.surface_names = sorted(all_surface_names)

        idx = self.surface_selector.FindString(selected_name)
        if idx == wx.NOT_FOUND:
            idx = 0

        self.surface_selector.Set(self.surface_names)
        self.surface_selector.SetSelection(idx)

        self.active_surface = None
        self.need_update = True

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

        # Choose surface to plot
        h = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h, flag=wx.ALL | wx.EXPAND, border=4)
        h.Add(wx.StaticText(self, label='Surface:'), flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(5)
        self.surface_selector = wx.Choice(self)
        self.surface_selector.Bind(wx.EVT_CHOICE, self.reset_surfaces)
        h.Add(self.surface_selector, proportion=1)
        h.AddSpacer(5)
        b = wx.Button(self, label='Load')
        b.Bind(wx.EVT_BUTTON, self.load_surfaces)
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
        fgs.Add(
            wx.StaticText(self, label='Aspect ratio:'), flag=wx.ALIGN_CENTER_VERTICAL
        )
        self.aspect = wx.Choice(self)
        self.aspect.Set(['equal', 'auto'])
        self.aspect.Select(0)
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

    def load_surfaces(self, evt=None):
        isurf = self.surface_selector.GetSelection()
        surface_name = self.surface_names[isurf]
        self.active_surface = surface_name

        max_num_timesteps = 0
        tmin, tmax = 1e100, -1e100
        xmin, ymin, xmax, ymax = 1e100, 1e100, -1e100, -1e100

        # Populate the cache for these surfaces
        for results in self.istate.active_results:
            if surface_name in results.surfaces:
                surf = results.surfaces[surface_name]

                with wx.BusyCursor():
                    try:
                        _description, _value, _dim, timesteps, data = surf.get_surfaces(
                            cache=True
                        )
                    except BaseException:
                        continue

                    # Compute bounds
                    tmin = min(tmin, timesteps[0])
                    tmax = max(tmax, timesteps[-1])
                    max_num_timesteps = max(max_num_timesteps, len(timesteps))
                    for contours in data:
                        for contour in contours:
                            xmin = min(xmin, numpy.min(contour[0]))
                            ymin = min(ymin, numpy.min(contour[1]))
                            xmax = max(xmax, numpy.max(contour[0]))
                            ymax = max(ymax, numpy.max(contour[1]))

        self.tmin = tmin
        self.tmax = tmax
        self.max_num_timesteps = max_num_timesteps
        self.time_selector.SetMax(self.max_num_timesteps)

        xdiff = xmax - xmin
        ydiff = ymax - ymin
        xmin, xmax = xmin - 0.1 * xdiff, xmax + 0.1 * xdiff
        ymin, ymax = ymin - 0.1 * ydiff, ymax + 0.1 * ydiff
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.title.ChangeValue('Iso surface %s (t=%%g)' % surface_name)
        self.xlabel.ChangeValue('x')
        self.ylabel.ChangeValue('y')

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
        if self.active_surface is None:
            self.canvas.draw()
            return
        surface_name = self.active_surface

        # Get time
        f = self.time_selector.GetValue() * 1.0 / self.max_num_timesteps
        t = self.tmin + (self.tmax - self.tmin) * f

        tsel = set()
        for ires, results in enumerate(self.istate.active_results):
            if surface_name in results.surfaces:
                surf = results.surfaces[surface_name]
                with wx.BusyCursor():
                    try:
                        _description, _value, _dim, timesteps, data = surf.get_surfaces(
                            cache=True
                        )
                    except Exception:
                        continue

                i = numpy.argmin(abs(timesteps - t))
                dt = timesteps[1] - timesteps[0]
                if abs(timesteps[i] - t) > 1.5 * dt:
                    continue
                tsel.add(timesteps[i])

                xvals = []
                yvals = []
                for contour in data[i]:
                    if xvals:
                        xvals.append(None)
                        yvals.append(None)
                    xvals.extend(contour[0])
                    yvals.extend(contour[1])
                self.axes.plot(xvals, yvals, color='C%d' % ires, label=results.label)

        self.axes.set_xlim(self.xmin, self.xmax)
        self.axes.set_ylim(self.ymin, self.ymax)
        self.axes.set_aspect(self.aspect.GetStringSelection())

        if len(tsel) == 1:
            t = tsel.pop()

        title = self.title.GetValue()
        try:
            title = title % t
        except Exception:
            pass

        self.axes.set_title(title)
        self.axes.set_xlabel(self.xlabel.GetValue())
        self.axes.set_ylabel(self.ylabel.GetValue())
        if len(self.istate.results) > 1:
            self.axes.legend(loc='lower right')
        self.fig.tight_layout()

        self.canvas.draw()
        self.need_update = False
