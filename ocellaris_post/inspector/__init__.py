try:
    import wx
except ImportError:
    print 'Missing python module "wx"'
    print
    print 'You must install wxPython to run the GUI. Python wheels are'
    print 'available for most platforms in addition to conda and other'
    print 'packages. The code has been tested with wxPython-4.0.0a3, an'
    print 'alpha release of wxPython 4 (which seems to work perfectly).'
    print
    print 'ERROR: missing wxPython'
    exit(1)

from ocellaris_post import Results
from wx.lib.pubsub import pub


# PubSub topics
TOPIC_METADATA = 'updated_metadata'
TOPIC_RELOAD = 'reloaded_data'


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
