import wx
from ocellaris_post import Results
from wx.lib.pubsub import pub


# PubSub topics
TOPIC_METADATA = 'metadata_updated'


from .inspector import OcellarisInspector


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
