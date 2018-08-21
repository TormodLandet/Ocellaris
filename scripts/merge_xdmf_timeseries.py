import os
import sys


def merge_xdmf_timeseries(inp_file_names, out_file_name, verbose=True):
    """
    Merge multiple XDMF time series files

    Each input XDMF file can contain multiple time steps, which should be
    unique among all the input files. The time steps are sorted and
    written to a new XDMF file

    The current code is quite hacky and probably only works for input XDMF
    files written by the current (2018.1) version of DOLFIN
    """
    timesteps = []
    for fn in inp_file_names:
        if verbose:
            print('Reading %s' % fn)

        state = ''
        first_grid = None
        for line in open(fn, 'rt'):
            if state == '':
                if 'CollectionType="Temporal"' in line:
                    state = 'read_timesteps'
            elif state == 'read_timesteps':
                if '<Grid' in line:
                    state = 'read_timestep'
                    ts_data = line
                if '</Grid>' in line:
                    state = ''
            elif state == 'read_timestep':
                ts_data += line
                if '<Time' in line and 'Value=' in line:
                    t = line.split('Value="')[1].split('"')[0]
                    t = float(t)
                if '</Grid>' in line:
                    timesteps.append((t, first_grid, ts_data))
                    if first_grid is None:
                        first_grid = ts_data
                    state = 'read_timesteps'

    timesteps.sort()

    # Correct xi:include indices
    timesteps2 = []
    indices = {}
    i = 0
    for _, first_grid, grid in timesteps:
        i += 1
        if first_grid is None:
            indices[grid] = i
        else:
            idx = indices[first_grid]
            grid = grid.replace('Grid[1]', 'Grid[%d]' % idx)
        timesteps2.append(grid)

    if verbose:
        print(
            '\nFound %d time steps (from t=%r to t=%r)'
            % (len(timesteps), timesteps[0][0], timesteps[-1][0])
        )
        print('Writing %s' % out_file_name)

    with open(out_file_name, 'wt') as out:
        out.write('<?xml version="1.0"?>\n')
        out.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        out.write('<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        out.write('  <Domain>\n')
        out.write('    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
        for ts_data in timesteps2:
            out.write(ts_data)
        out.write('    </Grid>\n')
        out.write('  </Domain>\n')
        out.write('</Xdmf>\n')


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
            error("Input file %r does not exist!" % fn)

    if os.path.exists(out_file_name):
        error("Output file %r exists!" % out_file_name)

    merge_xdmf_timeseries(inp_file_names, out_file_name)
