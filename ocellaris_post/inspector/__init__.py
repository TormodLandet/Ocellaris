try:
    import wx
except ImportError:
    print('Missing python module "wx"')
    print()
    print('You must install wxPython to run the GUI. Python wheels are')
    print('available for most platforms in addition to conda and other')
    print('packages. The code has been tested with wxPython-4.0.0a3, an')
    print('alpha release of wxPython 4 (which seems to work perfectly).')
    print()
    print('ERROR: missing wxPython')
    exit(1)


import os
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


class InspectorState(object):
    def __init__(self):
        """
        Store the data to be inspected
        """
        self.results = []
        self.persistence = InspectorPersistence(self)
        
    @property
    def active_results(self):
        return [r for r in self.results if r.active_in_gui]
    
    def open(self, file_name, label=None):
        """
        Open a new result file
        """
        r = Results(file_name)
        self.persistence.set_label(r, label)
        self.results.append(r)
        r.active_in_gui = True
    
    def reload(self, only_active=True):
        """
        Reload the data. Usefull when plotting log files that are 
        continuously updated
        """
        for r in self.results:
            if r.active_in_gui or not only_active: 
                r.reload()
    
    def close(self, idx):
        """
        Close the results file with the given index
        """
        del self.results[idx]


class InspectorPersistence(object):
    def __init__(self, inspector_state):
        """
        Store some data between runs of the inspector so that the
        user does not have to reconfigure the program each time it 
        is started
        """
        self.istate = inspector_state
        
        # Cache dir per the "XDG Base Directory Specification"
        cache_dir_default = os.path.expanduser('~/.cache')
        cache_dir = os.environ.get('XDG_CACHE_HOME', cache_dir_default)
        self.cache_file_name = os.path.join(cache_dir, 'ocellaris_post_inspector.yaml')
        
        # Automatic saving a while after each metadata update
        pub.subscribe(self.save_soon, TOPIC_METADATA)
        self.timer = None
        
        self.load()
    
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
    
    def set_label(self, res, label):
        # Use label if provided
        if label is not None:
            res.label = label
            return
        
        # Get persisent label if it exists or default to file name as label
        lables = self._cached_data.setdefault('result_file_lables', {})        
        if res.file_name in lables:
            res.label = lables[res.file_name]
        else:
            res.label = os.path.basename(res.file_name)
    
    def save(self, evt=None):
        # Save lables
        lables = self._cached_data.setdefault('result_file_lables', {})
        for res in self.istate.results:
            assert res.label is not None
            lables[res.file_name] = res.label
        
        with open(self.cache_file_name, 'wb') as f:
            yaml.safe_dump(self._cached_data, f)
        
        self.timer = None


def setup_yaml():
    """
    Make PyYaml load and store keys in dictionaries 
    ordered like they were on the input file
    """
    mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

    def dict_representer(dumper, data):
        return dumper.represent_dict(data.items())

    def dict_constructor(loader, node):
        return collections.OrderedDict(loader.construct_pairs(node))

    yaml.add_representer(collections.OrderedDict, dict_representer)
    yaml.add_constructor(mapping_tag, dict_constructor)


def show_inspector(file_names, lables):
    """
    Show wxPython window that allows chosing which report to show
    """
    setup_yaml()
    
    istate = InspectorState()
    for file_name, label in zip(file_names, lables):
        if not os.path.isfile(file_name):
            raise IOError('The results file %r does not exist' % file_name)
        istate.open(file_name, label)
    
    app = wx.App()
    frame = OcellarisInspector(istate)
    frame.Show()
    app.MainLoop()
