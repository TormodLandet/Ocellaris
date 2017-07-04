import wx
from . import pub, TOPIC_METADATA, TOPIC_RELOAD, TOPIC_NEW_ACCEL


class OcellarisSetupPanel(wx.Panel):
    def __init__(self, parent, results):
        super(OcellarisSetupPanel, self).__init__(parent)
        self.results = results
        self.layout_widgets()
    
    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)
        
        #######################################################################
        st = wx.StaticText(self, label='Ocellaris simulation result lables:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        
        # Metadata FlexGridSizer
        Nrows = len(self.results)
        fgs = wx.FlexGridSizer(rows=Nrows, cols=3, vgap=3, hgap=10)
        fgs.AddGrowableCol(1, proportion=1)
        v.Add(fgs, flag=wx.ALL|wx.EXPAND, border=6)
        
        # Customize the lables
        self.label_controls = []
        for il, results in enumerate(self.results):
            st = wx.StaticText(self, label='File %d:' % il)
            st.SetToolTip(results.file_name)
            fgs.Add(st, flag=wx.ALIGN_CENTER_VERTICAL)
            
            label_ctrl = wx.TextCtrl(self, value=results.label)
            label_ctrl.Bind(wx.EVT_TEXT, self.update_lables)
            fgs.Add(label_ctrl, flag=wx.EXPAND)
            self.label_controls.append(label_ctrl)
            
            st = wx.StaticText(self, label='(%s)' % results.label)
            st.SetToolTip(results.file_name)
            fgs.Add(st, flag=wx.ALIGN_CENTER_VERTICAL)
        
        #######################################################################
        st = wx.StaticText(self, label='Reload running simulation data:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        b = wx.Button(self, label='Reload (Ctrl+R)')
        b.Bind(wx.EVT_BUTTON, self.reload_data)
        v.Add(b, flag=wx.ALL, border=10)
        pub.sendMessage(TOPIC_NEW_ACCEL, callback=self.reload_data, key='R')
        
        v.Fit(self)
    
    def update_lables(self, event):
        for label_control, results in zip(self.label_controls, self.results):
            results.label = label_control.GetValue()
        pub.sendMessage(TOPIC_METADATA)
    
    def reload_data(self, evt=None):
        with wx.BusyCursor():
            for results in self.results:
                results.reload()
            pub.sendMessage(TOPIC_METADATA)
            pub.sendMessage(TOPIC_RELOAD)
