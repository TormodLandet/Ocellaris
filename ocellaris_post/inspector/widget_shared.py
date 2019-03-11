# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import wx

# For code eval in PlotCustomLine
from math import *
from numpy import *


class PlotLimSelectors(wx.Panel):
    def __init__(self, parent, callback):
        super().__init__(parent)
        self.callback = callback
        self.layout_widgets()
        self.timer = None

    def layout_widgets(self):
        fgs = wx.FlexGridSizer(rows=1, cols=10, vgap=3, hgap=9)
        for i in range(4):
            fgs.AddGrowableCol(1 + i * 2, proportion=1)
        self.SetSizer(fgs)

        self.ctrls = []
        for label in 'x-min x-max y-min y-max'.split():
            fgs.Add(wx.StaticText(self, label=label), flag=wx.ALIGN_CENTER_VERTICAL)
            ctrl = wx.TextCtrl(self, size=(20, -1))
            ctrl.Bind(wx.EVT_TEXT, self.callback_soon)
            fgs.Add(ctrl, flag=wx.EXPAND)
            self.ctrls.append(ctrl)

        b = wx.Button(self, label='Clear')
        b.Bind(wx.EVT_BUTTON, self.clear_lim)
        fgs.Add(b, flag=wx.EXPAND)

    def clear_lim(self, evt=None):
        for ctrl in self.ctrls:
            ctrl.Value = ''

    def callback_soon(self, evt=None):
        if self.timer is None:
            self.timer = wx.CallLater(500, self.callback_now)

    def callback_now(self, evt=None):
        self.timer = None
        self.callback()

    def _get_lim(self, c0, c1):
        try:
            minval = float(self.ctrls[c0].Value)
        except ValueError:
            minval = None

        try:
            maxval = float(self.ctrls[c1].Value)
        except ValueError:
            maxval = None

        return minval, maxval

    def get_xlim(self):
        return self._get_lim(0, 1)

    def get_ylim(self):
        return self._get_lim(2, 3)


class PlotCustomLine(wx.Panel):
    def __init__(self, parent, callback):
        super().__init__(parent)
        self.callback = callback
        self.layout_widgets()
        self.timer = None

    @property
    def active(self):
        return self.cb_active.Value

    def layout_widgets(self):
        fgs = wx.FlexGridSizer(rows=1, cols=5, vgap=3, hgap=9)
        for i in range(2):
            fgs.AddGrowableCol(1 + i * 2, proportion=2)
        self.SetSizer(fgs)

        self.ctrls = []
        for label, value in (
            ('Function', '2*x**2'),
            ('Name', 'Custom function y=2*x^2'),
        ):
            fgs.Add(wx.StaticText(self, label=label), flag=wx.ALIGN_CENTER_VERTICAL)
            ctrl = wx.TextCtrl(self, size=(20, -1), value=value)
            ctrl.Bind(wx.EVT_TEXT, self.callback_soon)
            fgs.Add(ctrl, flag=wx.EXPAND)
            self.ctrls.append(ctrl)

        self.cb_active = wx.CheckBox(self, label='Active')
        self.cb_active.Value = False
        self.cb_active.Bind(wx.EVT_CHECKBOX, self.callback_soon)
        fgs.Add(self.cb_active, flag=wx.EXPAND)

    def callback_soon(self, evt=None):
        if self.timer is None:
            self.timer = wx.CallLater(500, self.callback_now)

    def callback_now(self, evt=None):
        self.timer = None
        self.callback()

    def get_function(self, x):
        N = len(x)
        code = self.ctrls[0].Value
        name = self.ctrls[1].Value
        try:
            y = eval(code)
            assert len(y) == N
            return y, name
        except Exception:
            return x * 0 - 1, name
