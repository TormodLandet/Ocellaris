# encoding: utf8
import os
import collections
import re
import yaml


def read_yaml_input_file(file_name=None, yaml_string=None, error_func=None):
    """
    Read the input to an Ocellaris simulation from a YAML formated input file or a
    YAML formated string. The user will get an error if the input is malformed

    This reader is less advanced than what is used inside Ocellaris, but it does
    not depend on FEniCS or any heavy libraries, so it is more lightweight. For
    full compatibility the ``ocellaris.Simulation().input.read_yaml(...)`` method
    should be used.
    """
    if error_func is None:

        def error_func(header, description):
            raise ValueError('%s\n%s' % (header, description))

    setup_yaml()

    def load_ocellaris_input(file_name=None, yaml_string=None):
        if yaml_string is None:
            with open(file_name, 'rt') as inpf:
                yaml_string = inpf.read()
        else:
            assert file_name is None

        try:
            inp = yaml.load(yaml_string)
        except ValueError as e:
            error_func('Error on input file', str(e))
        except yaml.YAMLError as e:
            error_func('Input file "%s" is not a valid YAML file' % file_name, str(e))

        assert 'ocellaris' in inp
        assert inp['ocellaris']['type'] == 'input'
        assert inp['ocellaris']['version'] == 1.0
        return inp

    inp_dicts = []
    seen_files = set()

    def add_input_dictionary_to_list(inp):
        inp_dicts.append(inp)

        for base_file_name in inp['ocellaris'].get('bases', []):
            fn = search_for_file(base_file_name, file_name, error_func)

            # Make sure there are no circular dependencies
            fn = os.path.abspath(fn)
            assert fn not in seen_files, 'Circular dependency in input file bases'
            seen_files.add(fn)

            # Load base input file
            base_inp = load_ocellaris_input(fn)
            add_input_dictionary_to_list(base_inp)

    # Read the YAML data, potentially with base input files
    inp = load_ocellaris_input(file_name, yaml_string)
    add_input_dictionary_to_list(inp)
    final = collections.OrderedDict()
    for inp in inp_dicts[::-1]:
        merge_nested_dicts(final, inp)
    return final


def setup_yaml():
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
    yaml.reader.Reader.NON_PRINTABLE = re.compile(
        '[^\x09\x0A\x0D\x20-\x7E\x85\xA0-\uD7FF\uE000-\uFFFD\U00010000-\U0010FFFF]'
    )


def search_for_file(file_name, inp_file_name, error_func):
    """
    Serch first relative to the current working dir and then
    relative to the input file dir
    """
    # Check if the path is absolute or relative to the
    # working directory
    if os.path.exists(file_name):
        return file_name

    # Check if the path is relative to the input file dir
    inp_file_dir = os.path.dirname(inp_file_name)
    pth2 = os.path.join(inp_file_dir, file_name)
    if os.path.exists(pth2):
        return pth2

    error_func('File not found', 'The specified file "%s" was not found' % file_name)


def merge_nested_dicts(base, child):
    """
    Merge child dictionary into parent dictionary, overwriting and adding data
    from the child, keeping the base value if the child value is not specified
    for a given key. Dictionaries are merged recursively, other data like lists
    are treated as values and just replaced directly if found.
    """
    for k, v in child.items():
        if (
            not k in base
            or not isinstance(base[k], collections.abc.Mapping)
            or not isinstance(v, collections.abc.Mapping)
        ):
            base[k] = v
        else:
            merge_nested_dicts(base[k], v)
