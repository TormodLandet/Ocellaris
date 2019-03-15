# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os


MY_DIR = os.path.abspath(os.path.dirname(__file__))
SRC_DIR = os.path.join(MY_DIR, '..')
SCHEMA = 'ocellaris/input_file_schema.yml'


def pytest_generate_tests(metafunc):
    """
    Setup input files to test (all *.inp files in the same
    directory as this file
    """
    def find_inp_files(dirname):
        path = os.path.join(SRC_DIR, dirname)
        for fn in os.listdir(path):
            pth_rel = os.path.join(dirname, fn)
            pth_abs = os.path.join(path, fn)
            if os.path.isdir(pth_abs):
                yield from find_inp_files(pth_rel)
            if fn.endswith('inp'):
                yield pth_rel[2:]

    inpfiles = list(find_inp_files('./demos'))
    metafunc.parametrize("inp_file", inpfiles)


def test_valid_input_file(inp_file):
    """
    Make sure all input files in the Ocellaris repository can be read by
    the input file reader and pass YAML schema validation by yschema
    """
    import collections
    import yaml
    import yschema

    # Read input files in an ordered way
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         lambda loader, node: collections.OrderedDict(loader.construct_pairs(node)))

    data_file = os.path.join(SRC_DIR, inp_file)
    schema_file = os.path.join(SRC_DIR, SCHEMA)

    with open(data_file, 'rt') as f:
        data = yaml.load(f)

    with open(schema_file, 'rt') as f:
        schema = yaml.load(f)

    yschema.validate(data, schema)
