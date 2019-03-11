# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from .results import Results
from .input import read_yaml_input_file  # noqa


def open_results(file_or_dir, derived=True, inner_iterations=True) -> Results:
    """
    Give a file or directory name and get back a Results object. File names
    must point to an Ocellaris log file or an Ocellaris restart file.
    Directory names can be given for directories that contain only one log or
    restart file.

    If derived is True then derived data, such as the timestep dt and the total
    energy Et will be included. Otherwise only data explicityly written to an
    output file will be included.

    If inner_iterations is True then data from sub-timestep iterations will also
    be included, otherwise only data at the end of the time step is included.
    """
    import os
    from .files import get_result_file_name

    if os.path.isdir(file_or_dir):
        file_or_dir = get_result_file_name(file_or_dir)

    if os.path.isfile(file_or_dir):
        return Results(file_or_dir, derived=derived, inner_iterations=inner_iterations)
    else:
        print('ERROR: not a file %r' % file_or_dir)
