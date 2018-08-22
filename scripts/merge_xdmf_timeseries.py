import os
import sys
import h5py
from xml.etree import ElementTree as ET


def merge_xdmf_timeseries(inp_file_names, out_file_name, verbose=True):
    """
    Merge multiple XDMF time series files

    Each input XDMF file can contain multiple time steps, which should be
    unique among all the input files. The time steps are sorted and
    written to a new XDMF file

    The current code is quite hacky and probably only works for input XDMF
    files written by the current (2018.1) version of DOLFIN
    """
    # Handle XML name spaces
    ns = {'xi': 'http://www.w3.org/2001/XInclude'}
    for k, v in ns.items():
        ET.register_namespace(k, v)

    # Collect grids (time steps with functions and mesh info)
    # The time series is also a grid, with GridType = "Collection" and
    # CollectionType = "Temporal", the time steps are sub-grids of this
    grids = []
    out_tree = None
    for fn in inp_file_names:
        if verbose:
            print('Reading %s' % fn)
        xml = ET.parse(fn)
        xdmf = xml.getroot()

        # Use the first file as a template for the output
        if out_tree is None:
            out_tree = ET.parse(fn)
            out_domain = out_tree.getroot().find('Domain')
            out_domain.clear()

        # Get the time step for this grid and store which previous
        # grid contains the mesh data (assumed to be the first grid)
        for tselem in xdmf.findall('.//Grid[@CollectionType="Temporal"]'):
            grid_with_mesh = None
            for gridelem in tselem.findall('Grid'):
                t = float(gridelem.find('Time[@Value]').attrib['Value'])
                grids.append((t, grid_with_mesh, gridelem))
                if grid_with_mesh is None:
                    grid_with_mesh = gridelem

    # Sort grids by simulation time step
    grids.sort()

    # Correct xi:include references
    grids2 = []
    i = 0
    indices = {}
    for t, grid_with_mesh, grid in grids:
        i += 1
        if grid_with_mesh is None:
            indices[grid] = i
        else:
            idx = indices[grid_with_mesh]
            inc = grid.find('xi:include', ns)
            ptr = inc.get('xpointer')
            inc.set('xpointer', ptr.replace('Grid[1]', 'Grid[%d]' % idx))
        grids2.append(grid)

    # Merge HDF5 files
    print('\nCopying datasets to merged HDF5 file')
    current_h5_file = source = None
    h5_out_file = out_file_name.replace('.xdmf', '.h5')
    with h5py.File(h5_out_file, 'w') as dest:
        counters = {}
        for grid in grids2:
            for ds in grid.findall('.//DataItem[@Format="HDF"]'):
                filename, pth = ds.text.rsplit(':', 1)
                prefix = 'meshdata' if 'Mesh' in pth else 'vec'
                counters[prefix] = counters.get(prefix, -1) + 1
                new_pth = '/%s_%05d' % (prefix, counters[prefix])

                if filename != current_h5_file:
                    if source is not None:
                        source.close()
                    print('Reading %s' % filename)
                    source = h5py.File(filename, 'r')
                    current_h5_file = filename

                source.copy(pth, dest, new_pth)
                ds.text = '%s:%s' % (h5_out_file, new_pth)
        source.close()
        print('Writing %s' % h5_out_file)

    if verbose:
        print('\nFound %d time steps (from t=%r to t=%r)' % (len(grids), grids[0][0], grids[-1][0]))
        print('Writing %s' % out_file_name)

    out_ts_elem = ET.SubElement(out_domain, 'Grid')
    out_ts_elem.set('Name', 'TimeSeries')
    out_ts_elem.set('GridType', 'Collection')
    out_ts_elem.set('CollectionType', 'Temporal')
    for grid in grids2:
        out_ts_elem.append(grid)
    out_tree.write(out_file_name, xml_declaration=True)


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
