import wx
import os
from . import pub, TOPIC_METADATA, TOPIC_RELOAD, TOPIC_NEW_ACCEL
from .dialog_cluster_connector import OcellarisClusterConnectorDialog

class OcellarisSetupPanel(wx.Panel):
    def __init__(self, parent, inspector_state):
        super(OcellarisSetupPanel, self).__init__(parent)
        self.istate = inspector_state
        self.layout_widgets()
    
    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)
        
        #######################################################################
        st = wx.StaticText(self, label='Open files:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        h = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h, flag=wx.EXPAND|wx.ALL, border=10)
        
        b = wx.Button(self, label='Open new file (Ctrl+O)')
        b.Bind(wx.EVT_BUTTON, self.select_file_to_open)
        pub.sendMessage(TOPIC_NEW_ACCEL, callback=self.select_file_to_open, key='O')
        h.Add(b)
        h.AddSpacer(5)
        b = wx.Button(self, label='Open running simulations')
        b.Bind(wx.EVT_BUTTON, self.show_cluster_connector)
        h.Add(b)
        
        #######################################################################
        st = wx.StaticText(self, label='Change lables of open files:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        
        st = wx.StaticText(self, label='Hover over "File X" label to see full file name')
        v.Add(st, flag=wx.ALL, border=8)
        
        # Metadata FlexGridSizer
        self.file_lable_sizer = wx.FlexGridSizer(rows=1, cols=4, vgap=3, hgap=10)
        self.file_lable_sizer.AddGrowableCol(1, proportion=1)
        v.Add(self.file_lable_sizer, flag=wx.ALL|wx.EXPAND, border=6)
        
        #######################################################################
        st = wx.StaticText(self, label='Reload running simulation data:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        b = wx.Button(self, label='Reload (Ctrl+R)')
        b.Bind(wx.EVT_BUTTON, self.reload_data)
        v.Add(b, flag=wx.ALL, border=10)
        pub.sendMessage(TOPIC_NEW_ACCEL, callback=self.reload_data, key='R')
                
        #######################################################################
        v.AddStretchSpacer(prop=1)
        h = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h, flag=wx.EXPAND|wx.ALL, border=10)
        h.Add(wx.StaticText(self, label='Open previous file:'), flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(10)
        self.prevfiles = wx.Choice(self)
        h.Add(self.prevfiles, proportion=1, flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(10)
        b = wx.Button(self, label='Open')
        b.Bind(wx.EVT_BUTTON, self.open_previous_file)
        h.Add(b, flag=wx.ALIGN_CENTER_VERTICAL)
        
        self.update_open_files()
        v.Fit(self)
    
    def update_open_files(self, _evt=None):
        fgs = self.file_lable_sizer
        fgs.Clear(delete_windows=True)
        fgs.SetRows(len(self.istate.results) + 1)
        
        # Customize the lables
        self.label_controls = []
        for il, results in enumerate(self.istate.results):
            st = wx.StaticText(self, label='File %d:' % il)
            st.SetToolTip(results.file_name)
            fgs.Add(st, flag=wx.ALIGN_CENTER_VERTICAL)
            
            label_ctrl = wx.TextCtrl(self, value=results.label)
            label_ctrl.Bind(wx.EVT_TEXT, self.update_lables)
            fgs.Add(label_ctrl, flag=wx.EXPAND)
            self.label_controls.append(label_ctrl)
            
            b = wx.Button(self, label="Close")
            b.SetToolTip('Close the results file %r' % results.file_name)
            
            def make_closer(il):
                def close(_evt=None):
                    with wx.BusyCursor():
                        self.istate.close(il)
                        pub.sendMessage(TOPIC_METADATA)
                        pub.sendMessage(TOPIC_RELOAD)
                        self.update_open_files()
                return close
            
            b.Bind(wx.EVT_BUTTON, make_closer(il))
            fgs.Add(b, flag=wx.ALIGN_CENTER_VERTICAL)
            
            def make_activator(il):
                def activate(evt):
                    with wx.BusyCursor():
                        self.istate.results[il].active_in_gui = evt.IsChecked() 
                        pub.sendMessage(TOPIC_METADATA)
                return activate
            
            tb = wx.ToggleButton(self, label='Active')
            tb.SetValue(results.active_in_gui)
            tb.Bind(wx.EVT_TOGGLEBUTTON, make_activator(il))
            fgs.Add(tb, flag=wx.ALIGN_CENTER_VERTICAL)
        
        self.GetSizer().Layout()
        
        # Update list of previous files
        self.prev_file_names = self.istate.persistence.get_prev_files(20)
        def shorten(fn, N=80):
            if len(fn) > N:
                fn = '...' + fn[3-N:]
            return fn
        self.prevfiles.Set([shorten(fn) for fn in self.prev_file_names])
        i = len(self.prev_file_names) - 1 - len(self.istate.results)
        if i >= 0:
            self.prevfiles.Select(i)
    
    def select_file_to_open(self, _evt=None):
        defdir = ''
        deffile = ''
        prev = self.istate.persistence.get_prev_files()
        if prev:
            defdir = os.path.split(prev[-1])[0] 
        
        wildcard = "H5 checkpoint and LOG files (*.h5;*.log)|*.h5;*.log"
        dlg = wx.FileDialog(self, "Open Occelaris results file", defdir, deffile,
                            wildcard, wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if dlg.ShowModal() != wx.ID_CANCEL:
            self.open_file(dlg.GetPath())
    
    def open_previous_file(self, _evt=None):
        idx = self.prevfiles.GetSelection()
        if idx == wx.NOT_FOUND:
            return
        file_name = self.prev_file_names[idx]
        self.open_file(file_name)
    
    def open_file(self, file_name):
        with wx.BusyCursor():
            self.istate.open(file_name)
            pub.sendMessage(TOPIC_METADATA)
            pub.sendMessage(TOPIC_RELOAD)
            self.update_open_files()
    
    def show_cluster_connector(self, _evt=None):
        with OcellarisClusterConnectorDialog(self, self.istate, self.open_file) as dlg:
            dlg.ShowModal()
    
    def update_lables(self, _evt=None):
        for label_control, results in zip(self.label_controls, self.istate.results):
            results.label = label_control.GetValue()
        pub.sendMessage(TOPIC_METADATA)
    
    def reload_data(self, _evt=None):
        with wx.BusyCursor():
            self.istate.reload()
            pub.sendMessage(TOPIC_METADATA)
            pub.sendMessage(TOPIC_RELOAD)
