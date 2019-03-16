# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
Inspect timestep reports from one or more Ocellaris restart files
"""
import wx
from . import pub, TOPIC_NEW_ACCEL
from .icons import OCELLARIS_ICON
from .panel_setup import OcellarisSetupPanel
from .panel_results import OcellarisReportsPanel
from .panel_files import OcellarisFilesPanel
from .panel_stairs import OcellarisStairsPanel
from .panel_surfaces import OcellarisSurfacesPanel
from .panel_line_probes import OcellarisLineProbePanel


class OcellarisInspector(wx.Frame):
    def __init__(self, inspector_state):
        super().__init__(None, title='Ocellaris Report Inspector')
        self.istate = inspector_state

        # Keyboard shortcuts
        self.accelerators = []
        pub.subscribe(self.add_accelerator, TOPIC_NEW_ACCEL)

        self.layout_widgets()
        self.SetSize(1000, 800)
        self.SetIcon(OCELLARIS_ICON.GetIcon())
        self.SetDropTarget(InspectorDropTarget(self, self.setup_panel.open_file))

    def layout_widgets(self):
        p = wx.Panel(self)
        nb = wx.Notebook(p, style=wx.NB_BOTTOM)

        self.setup_panel = OcellarisSetupPanel(nb, self.istate)
        nb.AddPage(self.setup_panel, 'Setup')
        self.setup_panel.SetBackgroundColour(p.GetBackgroundColour())

        self.files_panel = OcellarisFilesPanel(nb, self.istate)
        nb.AddPage(self.files_panel, 'Files')
        self.files_panel.SetBackgroundColour(p.GetBackgroundColour())

        self.reports_panel = OcellarisReportsPanel(nb, self.istate)
        nb.AddPage(self.reports_panel, 'Timestep reports')
        self.reports_panel.SetBackgroundColour(p.GetBackgroundColour())

        self.stairs_panel = OcellarisStairsPanel(nb, self.istate)
        nb.AddPage(self.stairs_panel, 'Stairs')
        self.stairs_panel.SetBackgroundColour(p.GetBackgroundColour())

        self.line_probes_panel = OcellarisLineProbePanel(nb, self.istate)
        nb.AddPage(self.line_probes_panel, 'Line probes')
        self.line_probes_panel.SetBackgroundColour(p.GetBackgroundColour())

        self.surfaces_panel = OcellarisSurfacesPanel(nb, self.istate)
        nb.AddPage(self.surfaces_panel, 'Surfaces')
        self.surfaces_panel.SetBackgroundColour(p.GetBackgroundColour())

        if self.istate.results:
            nb.SetSelection(2)
        else:
            nb.SetSelection(0)
        s = wx.BoxSizer()
        s.Add(nb, 1, wx.EXPAND)
        p.SetSizer(s)

    def add_accelerator(self, callback, key):
        """
        Add a keyboard shortcut, e.g. Ctrl+R for reload
        """
        new_id = wx.NewId()
        self.Bind(wx.EVT_MENU, callback, id=new_id)
        ae = wx.AcceleratorEntry()
        ae.Set(wx.ACCEL_CTRL, ord(key), new_id)
        self.accelerators.append(ae)

        atable = wx.AcceleratorTable(self.accelerators)
        self.SetAcceleratorTable(atable)


class InspectorDropTarget(wx.FileDropTarget):
    def __init__(self, inspector_frame, opener):
        """
        Support drag and drop of files to be visualised
        """
        super().__init__()
        self.opener = opener

    def OnDropFiles(self, x, y, file_names):
        for fn in file_names:
            self.opener(fn)

        return True  # We wanted the files
