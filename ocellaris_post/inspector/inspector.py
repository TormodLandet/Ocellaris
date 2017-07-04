"""
Inspect timestep reports from one or more Ocellaris restart files 
"""
import wx
from ocellaris_post import Results
from .icons import OCELLARIS_ICON
from .panel_results import OcellarisReportsPanel
from .panel_setup import OcellarisSetupPanel
from .panel_surfaces import OcellarisSurfacesPanel


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
        
        self.metadata_panel = OcellarisSetupPanel(nb, self.results)
        nb.AddPage(self.metadata_panel, 'Setup')
        self.metadata_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        self.reports_panel = OcellarisReportsPanel(nb, self.results)
        nb.AddPage(self.reports_panel, 'Timestep reports')
        self.reports_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        self.surfaces_panel = OcellarisSurfacesPanel(nb, self.results)
        nb.AddPage(self.surfaces_panel, 'Surfaces')
        self.surfaces_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        nb.SetSelection(1)
        s = wx.BoxSizer()
        s.Add(nb, 1, wx.EXPAND)
        p.SetSizer(s)


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
