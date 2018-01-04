import os, collections
from ocellaris.utils import ocellaris_error, get_root_value
import yaml


class UndefinedParameter(object):
    def __repr__(self):
        "For Sphinx"
        return '<UNDEFINED>'
UNDEFINED = UndefinedParameter()


class Input(collections.OrderedDict):
    def __init__(self, simulation, values=None, basepath=''):
        """
        Holds the input values provided by the user
        """
        if values:
            super(Input, self).__init__(values.items())
        else:
            super(Input, self).__init__()
        
        self.simulation = simulation
        self.basepath = basepath
        self._already_logged = set()
        
        if basepath and not basepath.endswith('/'):
            self.basepath = basepath + '/'
    
    def read_yaml(self, file_name=None, yaml_string=None):
        """
        Read the input to an Ocellaris simulation from a YAML formated input file or a 
        YAML formated string. The user will get an error if the input is malformed 
        """
        self._setup_yaml()
        
        if yaml_string is None:
            with open(file_name, 'rt') as inpf:
                yaml_string = inpf.read()
        else:
            assert file_name is None
        
        try:
            inp = yaml.load(yaml_string)
        except ValueError as e:
            ocellaris_error('Error on input file', str(e))
        except yaml.YAMLError as e:
            ocellaris_error('Input file "%s" is not a valid YAML file' % file_name, str(e))
        
        assert 'ocellaris' in inp
        assert inp['ocellaris']['type'] == 'input'
        assert inp['ocellaris']['version'] == 1.0
        
        self.clear()
        self.update(inp)
        self.file_name = file_name
    
    def get_value(self, path, default_value=UNDEFINED, required_type='any',
                  mpi_root_value=False, safe_mode=False):
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
                type="any" does no type checking. Other options
                are "int", "float", "string", "bool", "Input",
                "list(float)", "dict(string:any)" etc
            mpi_root_value: get the value on the root MPI process
            safe_mode: do not evaluate python expressions "py$ xxx"
        
        Returns:
            The input value if it exist otherwise the default value
        """
        # Allow path to be a list or a "/" separated string
        if isinstance(path, str):
            pathstr = self.basepath + path
            path = path.split('/')
        else:
            pathstr = self.basepath + '/'.join(path)
        
        def check_isinstance(value, classes):
            """
            Give error if the input data is not of the required type
            """
            value = eval_python_expression(self.simulation, value, pathstr, safe_mode)
            
            if not isinstance(value, classes):
                ocellaris_error('Malformed data on input file',
                                'Parameter %s should be of type %s,\nfound %r %r' % 
                                (pathstr, required_type, value, type(value)))
            return value
        
        def check_dict(d, keytype, valtype):
            """
            Check dict and eval any python expressions in the values
            """
            if d is None:
                # if every element in the dict is commented out the dict becomes None
                d = {}
            
            d = check_isinstance(d, dict_types)
            d_new = collections.OrderedDict()
            for key, val in d.items():
                check_isinstance(key, keytype)
                d_new[key] = check_isinstance(val, valtype)
            return d_new
        
        def check_list(d, valtype):
            """
            Check list and eval any python expressions in the values
            """
            d = check_isinstance(d, list)
            d_new = []
            for val in d:
                d_new.append(check_isinstance(val, valtype))
            return d_new
        
        # Get validation function according to required data type
        number = (int, int, float)
        dict_types = (dict, collections.OrderedDict)
        anytype = (int, float, str, list, tuple, dict, collections.OrderedDict, bool)
        if required_type == 'bool':
            def validate_and_convert(d): return check_isinstance(d, bool)
        elif required_type == 'float':
            # The YAML parser annoyingly thinks 1e-3 is a string (while 1.0e-3 is a float)
            def validate_and_convert(d):
                if isinstance(d, str):
                    try:
                        d = float(d)
                    except ValueError:
                        pass
                return check_isinstance(d, number)
        elif required_type == 'int':
            def validate_and_convert(d): return check_isinstance(d, int)
        elif required_type == 'string':
            def validate_and_convert(d):
                d = check_isinstance(d, str)
                # SWIG does not like Python 2 Unicode objects
                return str(d)
        elif required_type == 'Input':
            def validate_and_convert(d):
                d = check_isinstance(d, dict_types)
                return Input(self.simulation, d, basepath=pathstr)
        elif required_type == 'dict(string:any)':
            def validate_and_convert(d): return check_dict(d, str, anytype)
        elif required_type == 'dict(string:dict)':
            def validate_and_convert(d): return check_dict(d, str, dict_types)
        elif required_type == 'dict(string:list)':
            def validate_and_convert(d): return check_dict(d, str, list)
        elif required_type == 'dict(string:float)':
            def validate_and_convert(d): return check_dict(d, str, number)
        elif required_type == 'list(float)':
            def validate_and_convert(d): return check_list(d, number)
        elif required_type == 'list(int)':
            def validate_and_convert(d): return check_list(d, int)
        elif required_type == 'list(string)':
            def validate_and_convert(d): return check_list(d, str)
        elif required_type == 'list(dict)':
            def validate_and_convert(d): return check_list(d, dict_types)
        elif required_type == 'any':
            def validate_and_convert(d): return check_isinstance(d, anytype)
        else:
            raise ValueError('Unknown required_type %s' % required_type)
        
        # Look for the requested key
        d = self
        for p in path:
            if isinstance(d, list):
                # This is a list, assume the key "p" is an integer position
                try:
                    p = int(p)
                except ValueError:
                    ocellaris_error('List index not integer',
                                    'Not a valid list index:  %s' % p)
            elif d is None or p not in d:
                # This is an empty dict or a dict missing the key "p"
                if default_value is UNDEFINED:
                    ocellaris_error('Missing parameter on input file',
                                    'Missing required input parameter:\n  %s' % pathstr)
                else:
                    msg  = '    No value set for "%s", using default value %r' % (pathstr, default_value)
                    if not msg in self._already_logged:
                        self.simulation.log.debug(msg)
                        self._already_logged.add(msg)
                    if required_type == 'Input':
                        default_value = Input(self.simulation, default_value)
                    return default_value
            d = d[p]
        
        # Validate the input data and convert to the requested type
        d = validate_and_convert(d)
        
        # Get the value on the root process
        if mpi_root_value:
            d = get_root_value(d)
        
        # Show what input values we use
        msg  = '    Input value "%s" set to %r' % (pathstr, d)
        if not msg in self._already_logged:
            self.simulation.log.debug(msg)
            self._already_logged.add(msg)
        
        return d
    
    def set_value(self, path, value):
        """
        Set an input value by its path in the input dictionary
        
        Arguments:
            path: a list of path components or the "/" separated
                path to the variable in the input dictionary
            value: the value to set
        
        """
        # Allow path to be a list or a "/" separated string
        if isinstance(path, str):
            path = path.split('/')
        
        d = self
        for p in path[:-1]:
            if isinstance(d, list):
                try:
                    p = int(p)
                except ValueError:
                    ocellaris_error('List index not integer',
                                    'Not a valid list index:  %s' % p)
            elif p not in d:
                d[p] = {}
            d = d[p]
        d[path[-1]] = value
    
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
        
        ocellaris_error('File not found', 'The specified file "%s" was not found' % file_name)
    
    def _setup_yaml(self):
        """
        Make PyYaml load and store keys in dictionaries 
        ordered like they were on the input file
        """
        _mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
    
        def dict_representer(dumper, data):
            return dumper.represent_dict(data.items())
    
        def dict_constructor(loader, node):
            return collections.OrderedDict(loader.construct_pairs(node))
    
        yaml.add_representer(collections.OrderedDict, dict_representer)
        yaml.add_constructor(_mapping_tag, dict_constructor)
        
        # PyYAML bugfix to be able to read, e.g., 𝜏𝜀𝜁𝜃
        # See https://stackoverflow.com/a/44875714
        import re
        yaml.reader.Reader.NON_PRINTABLE = re.compile(
            '[^\x09\x0A\x0D\x20-\x7E\x85\xA0-\uD7FF\uE000-\uFFFD\U00010000-\U0010FFFF]')
    
    def __str__(self):
        inp = collections.OrderedDict(self.items())
        return yaml.dump(inp, indent=4)


def eval_python_expression(simulation, value, pathstr, safe_mode=False):
    """
    We run eval with the math functions and user constants available on string
    values that are prefixed with "py$" indicating that they are dynamic
    expressions and not just static strings
    """
    if not isinstance(value, str) or not value.startswith('py$'):
        return value
    
    if safe_mode:
        ocellaris_error('Cannot have Python expression here',
                        'Not allowed to have Python expression here:  %s' % pathstr)
     
    # remove "py$" prefix
    expr = value[3:]
    
    # Build dictionary of locals for evaluating the expression    
    eval_locals = {}
    
    import math
    for name in dir(math):
        if not name.startswith('_'):
            eval_locals[name] = getattr(math, name)
    
    global_inp = simulation.input
    user_constants = global_inp.get_value('user_code/constants', {}, 'dict(string:float)',
                                          safe_mode=True)
    for name, const_value in user_constants.items():
        eval_locals[name] = const_value
    
    eval_locals['simulation'] = simulation
    eval_locals['t'] = eval_locals['time'] = simulation.time
    eval_locals['it'] = eval_locals['timestep'] = simulation.timestep
    eval_locals['dt'] = simulation.dt
    eval_locals['ndim'] = simulation.ndim
    
    try:
        value = eval(expr, globals(), eval_locals)
    except Exception:
        simulation.log.error('Cannot evaluate python code for %s' % pathstr)
        simulation.log.error('Python code is %s' % expr)
        raise
    return value
