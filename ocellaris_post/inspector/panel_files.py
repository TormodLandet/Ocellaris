import yaml
import wx
from . import pub, TOPIC_METADATA, TOPIC_RELOAD


class OcellarisFilesPanel(wx.Panel):
    def __init__(self, parent, results):
        super(OcellarisFilesPanel, self).__init__(parent)
        self.results = results
        self.need_update = True
        
        self.layout_widgets()
        self.Bind(wx.EVT_IDLE, self.on_idle)
        pub.subscribe(self.reload_soon, TOPIC_RELOAD)
        pub.subscribe(self.reload_soon, TOPIC_METADATA)
        
    def reload_soon(self, evt=None):
        self.need_update = True
    
    def layout_widgets(self):
        v = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(v)
        
        # Choose results to plot
        h = wx.BoxSizer(wx.HORIZONTAL)
        v.Add(h, flag=wx.ALL|wx.EXPAND, border=5)
        h.Add(wx.StaticText(self, label='File:'), flag=wx.ALIGN_CENTER_VERTICAL)
        h.AddSpacer(5)
        labels = [res.label for res in self.results]
        self.results_selector = wx.Choice(self, choices=labels)
        self.results_selector.Select(0)
        self.results_selector.Bind(wx.EVT_CHOICE, self.switch_file)
        h.Add(self.results_selector, proportion=1)
        
        # Monospaced font
        font = wx.Font(10, wx.TELETYPE, wx.NORMAL, wx.NORMAL, False)
        
        st = wx.StaticText(self, label='Input:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        self.input_file = wx.TextCtrl(self, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP)
        self.input_file.SetFont(font)
        v.Add(self.input_file, flag=wx.ALL|wx.EXPAND, proportion=1, border=10)
        
        st = wx.StaticText(self, label='Log:')
        st.SetFont(st.GetFont().Bold())
        v.Add(st, flag=wx.ALL, border=5)
        self.log_file = wx.TextCtrl(self, style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_DONTWRAP)
        self.log_file.SetFont(font)
        v.Add(self.log_file, flag=wx.ALL|wx.EXPAND, proportion=2, border=10)
        
        v.Fit(self)
    
    def on_idle(self, evt=None):
        if self.need_update:
            self.switch_file()
    
    def switch_file(self, evt=None):
        """
        Switch which simulation is active
        """
        i = self.results_selector.GetSelection()
        results = self.results[i]
        
        inp = yaml.dump(results.input, default_flow_style=False, default_style='')
        log = results.log
        
        if len(log) > 40000:
            log = (log[:20000] + 
                   '\n\n...\n... MIDDLE SECTION SKIPPED ...\n...\n\n' +
                   log[-20000:])
        
        self.input_file.SetValue(inp)
        self.log_file.SetValue(log)
        
        self.need_update = False
