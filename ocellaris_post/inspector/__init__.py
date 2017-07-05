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


import yaml
import collections
from ocellaris_post import Results
from wx.lib.pubsub import pub


# PubSub topics
TOPIC_METADATA = 'updated_metadata'
TOPIC_RELOAD = 'reloaded_data'
TOPIC_NEW_ACCEL = 'new_keyboard_shortcut'


# Must import the inspector after the definition of TOPIC_*
from .inspector import OcellarisInspector


def setup_yaml():
    """
    Make PyYaml load and store keys in dictionaries 
    ordered like they were on the input file
    """
    mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

    def dict_representer(dumper, data):
        return dumper.represent_dict(data.iteritems())

    def dict_constructor(loader, node):
        return collections.OrderedDict(loader.construct_pairs(node))

    yaml.add_representer(collections.OrderedDict, dict_representer)
    yaml.add_constructor(mapping_tag, dict_constructor)


def show_inspector(file_names, lables):
    """
    Show wxPython window that allows chosing which report to show
    """
    setup_yaml()
    
    results = []
    for file_name, label in zip(file_names, lables):
        res = Results(file_name)
        res.label = label
        results.append(res)

    app = wx.App()
    frame = OcellarisInspector(results)
    frame.Show()
    app.MainLoop()
