import wx
from . import pub, TOPIC_METADATA, TOPIC_RELOAD, TOPIC_NEW_ACCEL


class OcellarisSetupPanel(wx.Panel):
    def __init__(self, parent, inspector_state):
        super(OcellarisSetupPanel, self).__init__(parent)
        self.istate = inspector_state
        self.layout_widgets()
    
    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)
        
        #######################################################################
        st = wx.StaticText(self, label='Ocellaris simulation result lables:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        
        st = wx.StaticText(self, label='Hover over "File X" label to see full file name')
        v.Add(st, flag=wx.ALL, border=8)
        
        # Metadata FlexGridSizer
        self.file_lable_sizer = wx.FlexGridSizer(rows=1, cols=4, vgap=3, hgap=10)
        self.file_lable_sizer.AddGrowableCol(1, proportion=1)
        v.Add(self.file_lable_sizer, flag=wx.ALL|wx.EXPAND, border=6)
        self.update_open_files()
        
        #######################################################################
        st = wx.StaticText(self, label='Reload running simulation data:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        b = wx.Button(self, label='Reload (Ctrl+R)')
        b.Bind(wx.EVT_BUTTON, self.reload_data)
        v.Add(b, flag=wx.ALL, border=10)
        pub.sendMessage(TOPIC_NEW_ACCEL, callback=self.reload_data, key='R')
        
        v.Fit(self)
    
    def update_open_files(self, event=None):
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
                def close(evt=None):
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
        
        # Open button
        fgs.AddSpacer(1)
        b = wx.Button(self, label='Open new file (Ctrl+O)')
        b.Bind(wx.EVT_BUTTON, self.select_file_to_open)
        fgs.Add(b)
        fgs.AddSpacer(1)
        pub.sendMessage(TOPIC_NEW_ACCEL, callback=self.select_file_to_open, key='O')
        
        self.GetSizer().Layout()
    
    def select_file_to_open(self, evt=None):
        wildcard = "H5 checkpoint and LOG files (*.h5;*.log)|*.h5;*.log"
        dlg = wx.FileDialog(self, "Open Occelaris results file", "", "", wildcard, wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if dlg.ShowModal() != wx.ID_CANCEL:
            self.open_file(dlg.GetPath())
    
    def open_file(self, file_name):
        with wx.BusyCursor():
            self.istate.open(file_name)
            pub.sendMessage(TOPIC_METADATA)
            pub.sendMessage(TOPIC_RELOAD)
            self.update_open_files()
    
    def update_lables(self, event=None):
        for label_control, results in zip(self.label_controls, self.istate.results):
            results.label = label_control.GetValue()
        pub.sendMessage(TOPIC_METADATA)
    
    def reload_data(self, evt=None):
        with wx.BusyCursor():
            self.istate.reload()
            pub.sendMessage(TOPIC_METADATA)
            pub.sendMessage(TOPIC_RELOAD)
