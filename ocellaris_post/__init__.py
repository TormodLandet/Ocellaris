from .results import Results
from .input import read_yaml_input_file  # noqa


def open_results(file_or_dir, derived=True, inner_iterations=True) -> Results:
    """
    Give a file or directory name and get back a
    Results object
    """
    import os
    from .files import get_result_file_name

    if os.path.isdir(file_or_dir):
        file_or_dir = get_result_file_name(file_or_dir)

    if os.path.isfile(file_or_dir):
        return Results(file_or_dir, derived=derived, inner_iterations=inner_iterations)
    else:
        print('ERROR: not a file %r' % file_or_dir)
