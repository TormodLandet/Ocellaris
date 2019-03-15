import os
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
    out_tree = None
    file_grids = []
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
            tmax_file = -1e100
            for gridelem in tselem.findall('Grid'):
                t = float(gridelem.find('Time[@Value]').attrib['Value'])
                data = (t, grid_with_mesh, gridelem)
                if grid_with_mesh is None:
                    print('    File t0 = %g' % t)
                    grid_with_mesh = gridelem
                    fgrids = []
                    file_grids.append((t, fgrids))
                fgrids.append(data)
                tmax_file = max(tmax_file, t)
            print('    File tN = %g' % tmax_file)

    # Make sure the later files overwrite the earlier if they have the
    # same time steps. The later will be what was used further on
    file_grids.sort(key=lambda data: data[0])
    timestep_grids = {}
    for _, grids in file_grids:
        for t, grid_with_mesh, gridelem in grids:
            timestep_grids[t] = (grid_with_mesh, gridelem)

    # Sort grids by simulation time step
    grids = []
    for t in sorted(timestep_grids):
        grid_with_mesh, gridelem = timestep_grids[t]
        grids.append((t, grid_with_mesh, gridelem))

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
    h5_out_file = out_file_name.replace('.xdmf', '.h5')
    current_h5_file = source = None
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


def delete_existing_file(out_file_name):
    """
    Remove any previous merged XDMF files
    """
    h5_out_file = out_file_name.replace('.xdmf', '.h5')
    if os.path.isfile(out_file_name):
        told = -1
        with open(out_file_name, 'rt') as f:
            for line in f:
                if '<Time Value=' in line:
                    try:
                        tiold = float(line.split('"')[1])
                    except ValueError:
                        continue
                    told = max(told, tiold)
        print('Deleting existing %s with max time = %r' % (out_file_name, told))
        os.unlink(out_file_name)
    if os.path.isfile(h5_out_file):
        print('Deleting existing %s' % h5_out_file)
        os.unlink(h5_out_file)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    parser.add_argument('output_file')
    parser.add_argument('--delete-old', action='store_true')
    args = parser.parse_args()

    for fn in args.input_files:
        if not os.path.isfile(fn):
            parser.print_help()
            print("\nERROR: Input file %r does not exist!" % fn)

    if os.path.exists(args.output_file):
        if args.delete_old:
            delete_existing_file(args.output_file)
        else:
            parser.print_help()
            print("\nERROR: Output file %r exists!" % args.output_file)

    merge_xdmf_timeseries(args.input_files, args.output_file)
