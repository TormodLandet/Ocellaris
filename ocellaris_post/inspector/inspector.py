"""
Inspect timestep reports from one or more Ocellaris restart files 
"""
import os
import yaml
import wx
from . import pub, TOPIC_NEW_ACCEL, TOPIC_METADATA
from .icons import OCELLARIS_ICON
from .panel_results import OcellarisReportsPanel
from .panel_setup import OcellarisSetupPanel
from .panel_surfaces import OcellarisSurfacesPanel
from .panel_files import OcellarisFilesPanel


class OcellarisInspector(wx.Frame):
    def __init__(self, results):
        super(OcellarisInspector, self).__init__(None, title='Ocellaris Report Inspector')
        self.results = results
        
        # Load any configuration available for these results (e.g lables)
        self.persistence = OcellarisInspectorPersistence(results)
        self.persistence.load()
        
        # File names as labels for unlabeled results
        for res in results:
            if res.label is None:
                res.label = os.path.basename(res.file_name)
        
        # Keyboard shortcuts
        self.accelerators = []
        pub.subscribe(self.add_accelerator, TOPIC_NEW_ACCEL)
        
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
        
        self.files_panel = OcellarisFilesPanel(nb, self.results)
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


class OcellarisInspectorPersistence(object):
    def __init__(self, results):
        """
        Store some data between runs of the inspector so that the
        user does not have to reconfigure the program each time it 
        is started
        """
        self.results = results
        
        # Cache dir per the "XDG Base Directory Specification"
        cache_dir_default = os.path.expanduser('~/.cache')
        cache_dir = os.environ.get('XDG_CACHE_HOME', cache_dir_default)
        self.cache_file_name = os.path.join(cache_dir, 'ocellaris_post_inspector.yaml')
        
        # Automatic saving a while after each metadata update
        pub.subscribe(self.save_soon, TOPIC_METADATA)
        self.timer = None
    
    def save_soon(self, evt=None):
        if self.timer is not None:
            # Allready going to save
            return
        
        # Save after 5 second of inactivity (to avoid slowdowns in case 
        # there are multiple updates in a row, which is likely)
        self.timer = wx.CallLater(5000, self.save)
    
    def load(self):
        if not os.path.isfile(self.cache_file_name):
            self._cached_data = {}
            return
        
        with open(self.cache_file_name, 'rb') as f:
            self._cached_data = yaml.safe_load(f)
        
        # Load any lables that are unset
        lables = self._cached_data.setdefault('result_file_lables', {})
        for res in self.results:
            if res.label is not None:
                continue
            if res.file_name not in lables:
                continue
            res.label = lables[res.file_name]
    
    def save(self, evt=None):
        # Save lables
        lables = self._cached_data.setdefault('result_file_lables', {})
        for res in self.results:
            assert res.label is not None
            lables[res.file_name] = res.label
        
        with open(self.cache_file_name, 'wb') as f:
            yaml.safe_dump(self._cached_data, f)
        
        self.timer = None
