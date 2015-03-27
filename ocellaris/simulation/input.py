import os, collections
import yaml
from ocellaris.utils import report_error


class UndefinedParameter(object):
    def __repr__(self):
        "For Sphinx"
        return '<UNDEFINED>'
UNDEFINED = UndefinedParameter()


class Input(collections.OrderedDict):
    def __init__(self, simulation):
        """
        Holds the input values provided by the user
        """
        super(Input, self).__init__()
        self.simulation = simulation
    
    def read_yaml(self, file_name):
        """
        Read the input to an Ocellaris simulation from a YAML 
        formated input file. The user will get an error if the
        input file is malformed 
        """
        self._setup_yaml()
        try:
            with open(file_name, 'rt') as inpf:
                inp = yaml.load(inpf)
        except ValueError as e:
            report_error('Error on input file', str(e))
        
        assert 'ocellaris' in inp
        assert inp['ocellaris']['type'] == 'input'
        assert inp['ocellaris']['version'] == 1.0
        
        self.clear()
        self.update(inp)
        self.file_name = file_name
    
    def get_value(self, path, default_value=UNDEFINED, required_type='any'):
        """
        Get an input value by its path in the input dictionary
        
        Gives an error if there is no default value supplied
        and the  input variable does not exist
        
        Arguments:
            path: a list of path components or the "/" separated
                path to the variable in the input dictionary
            default_value: the value to return if the path does
                not exist in the input dictionary
            required_type: expected type of the variable. Giving 
                type="any" does no type checking
        
        Returns:
            The input value if it exist otherwise the default value
        """
        # Allow path to be a list or a "/" separated string
        if isinstance(path, basestring):
            pathstr = path
            path = path.split('/')
        else:
            pathstr = '/'.join(path)
        
        d = self
        for p in path:
            if p not in d:
                if default_value is UNDEFINED:
                    report_error('Missing parameter on input file',
                                 'Missing required input parameter:\n  %s' % pathstr)
                else:
                    return default_value
            d = d[p]
        
        def check_isinstance(value, classes):
            """
            Give error if the input data is not of the required type
            """
            if not isinstance(value, classes):
                report_error('Malformed data on input file',
                             'Parameter %s should be of type %s,\nfound %r %r' % 
                             (pathstr, required_type, value, type(value)))
        
        # Validate according to required data type
        number = (int, long, float)
        dict_types = (dict, collections.OrderedDict)
        if required_type == 'bool':
            check_isinstance(d, bool)
        elif required_type == 'float':
            check_isinstance(d, number)
        elif required_type == 'int':
            check_isinstance(d, int)
        elif required_type == 'string':
            check_isinstance(d, basestring)
            # SWIG does not like Python 2 Unicode objects
            d = str(d)
        elif required_type == 'dict(string:any)':
            check_isinstance(d, dict_types)
            for key, val in d.items():
                check_isinstance(key, basestring)
        elif required_type == 'dict(string:dict)':
            check_isinstance(d, dict_types)
            for key, val in d.items():
                check_isinstance(key, basestring)
                check_isinstance(val, dict_types)
        elif required_type == 'dict(string:list)':
            check_isinstance(d, dict_types)
            for key, val in d.items():
                check_isinstance(key, basestring)
                check_isinstance(val, list)
        elif required_type == 'list(float)':
            check_isinstance(d, list)
            for elem in d:
                check_isinstance(elem, number)
        elif required_type == 'list(dict)':
            check_isinstance(d, list)
            for elem in d:
                check_isinstance(elem, dict_types)
        elif required_type == 'any':
            pass
        else:
            raise ValueError('Unknown required_type %s' % required_type)
        return d
    
    def get_output_file_path(self, path, default_value=UNDEFINED):
        """
        Get the name of an output file
        
        Automatically prefixes the file name with the output prefix
        """
        prefix = self.get_value('output/prefix', '')
        filename = self.get_value(path, default_value, 'string')
        if default_value is None and filename is None:
            return None
        else:
            return prefix + filename
        
    def get_input_file_path(self, file_name):
        """
        Serch first relative to the current working dir and then
        relative to the input file dir
        """
        # Check if the path is absolute or relative to the
        # working directory
        if os.path.exists(file_name):
            return file_name
        self.simulation.log.debug('File does not exist: %s' % file_name)
        
        # Check if the path is relative to the inouf file dir
        inp_file_dir = os.path.dirname(self.file_name)
        pth2 = os.path.join(inp_file_dir, file_name)
        if os.path.exists(pth2):
            return pth2
        self.simulation.log.debug('File does not exist: %s' % pth2)
        
        report_error('File not found', 'The specified file "%s" was not found' % file_name)
    
    def _setup_yaml(self):
        """
        Make PyYaml load and store keys in dictionaries 
        ordered like they were on the input file
        """
        _mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
    
        def dict_representer(dumper, data):
            return dumper.represent_dict(data.iteritems())
    
        def dict_constructor(loader, node):
            return collections.OrderedDict(loader.construct_pairs(node))
    
        yaml.add_representer(collections.OrderedDict, dict_representer)
        yaml.add_constructor(_mapping_tag, dict_constructor)
        
    def __str__(self):
        inp = collections.OrderedDict(self.items())
        return yaml.dump(inp, indent=4)