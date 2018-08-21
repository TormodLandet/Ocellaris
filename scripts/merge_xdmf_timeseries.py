import os, sys


def merge_xdmf_timeseries(inp_file_names, out_file_name):
    """
    Merge multiple XDMF time series files

    Each input XDMF file can contain multiple time steps, which should be
    unique among all the input files. The time steps are sorted and written
    to a new XDMF file
    """
    pass


def error(msg):
    print('USAGE:')
    print('    merge_xdmf_timeseries inp1.xdmf [inp2.xdmf ...] out.xdmf')
    print('ERROR!')
    print('   ', msg)
    sys.exit(1)


if __name__ == '__main__':
    fns = sys.argv[1:]
    assert len(fns) > 1, "You must give input file names (min 1) and then output file name"

    inp_file_names = fns[:-1]
    out_file_name = fns[-1]

    for fn in inp_file_names:
        if not os.path.isfile(fn):
            error("File name %r does not exist!" % fn)

    if os.path.exists(out_file_name):
        error("Output file exists!")

    merge_xdmf_timeseries(inp_file_names, out_file_name)
