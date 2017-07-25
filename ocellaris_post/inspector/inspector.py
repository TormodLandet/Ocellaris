"""
Inspect timestep reports from one or more Ocellaris restart files 
"""
import wx
from . import pub, TOPIC_NEW_ACCEL
from .icons import OCELLARIS_ICON
from .panel_results import OcellarisReportsPanel
from .panel_setup import OcellarisSetupPanel
from .panel_surfaces import OcellarisSurfacesPanel
from .panel_files import OcellarisFilesPanel


class OcellarisInspector(wx.Frame):
    def __init__(self, inspector_state):
        super(OcellarisInspector, self).__init__(None, title='Ocellaris Report Inspector')
        self.istate = inspector_state
        
        # Keyboard shortcuts
        self.accelerators = []
        pub.subscribe(self.add_accelerator, TOPIC_NEW_ACCEL)
        
        self.layout_widgets()
        self.SetSize(800, 800)
        self.SetIcon(OCELLARIS_ICON.GetIcon())
        self.SetDropTarget(InspectorDropTarget(self, self.setup_panel.open_file))
    
    def layout_widgets(self):
        p = wx.Panel(self)
        nb = wx.Notebook(p, style=wx.NB_BOTTOM)
        
        self.setup_panel = OcellarisSetupPanel(nb, self.istate)
        nb.AddPage(self.setup_panel, 'Setup')
        self.setup_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        self.reports_panel = OcellarisReportsPanel(nb, self.istate)
        nb.AddPage(self.reports_panel, 'Timestep reports')
        self.reports_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        self.surfaces_panel = OcellarisSurfacesPanel(nb, self.istate)
        nb.AddPage(self.surfaces_panel, 'Surfaces')
        self.surfaces_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        self.files_panel = OcellarisFilesPanel(nb, self.istate)
        nb.AddPage(self.files_panel, 'Files')
        self.files_panel.SetBackgroundColour(p.GetBackgroundColour())
        
        nb.SetSelection(1)
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
        super(InspectorDropTarget, self).__init__()
        self.opener = opener
    
    def OnDropFiles(self, x, y, file_names):
        for fn in file_names:
            self.opener(fn)
        
        return True # We wanted the files
