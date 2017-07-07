import wx


class PlotLimSelectors(wx.Panel):
    def __init__(self, parent, callback):
        super(PlotLimSelectors, self).__init__(parent)
        self.callback = callback
        self.layout_widgets()
        self.timer = None
    
    def layout_widgets(self):
        fgs = wx.FlexGridSizer(rows=1, cols=10, vgap=3, hgap=9)
        for i in range(4):
            fgs.AddGrowableCol(1 + i*2, proportion=1)
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
    
    def get_xlim(self):
        try: minval = float(self.ctrls[0].Value)
        except ValueError: minval = None
        
        try: maxval = float(self.ctrls[1].Value)
        except ValueError: maxval = None
        
        return minval, maxval
    
    def get_ylim(self):
        try: minval = float(self.ctrls[2].Value)
        except ValueError: minval = None
        
        try: maxval = float(self.ctrls[3].Value)
        except ValueError: maxval = None
        
        return minval, maxval
