# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import re
import yaml
import numpy
import h5py
import dolfin
from ocellaris.utils import ocellaris_error


class RestartFileIO:
    def __init__(self, simulation, persisted_python_data):
        """
        This class handles reading and writing the simulation state such as
        velocity and presure fields. Files for postprocessing (xdmf) are also
        handled here
        """
        self.simulation = simulation
        self.persisted_python_data = persisted_python_data

    def is_restart_file(self, file_name):
        """
        Is the given file an Ocellaris restart file
        """
        HDF5_SIGNATURE = b'\211HDF\r\n\032\n'
        try:
            # The HDF5 header is not guaranteed to be at offset 0, but for our
            # purposes this can be assumed as we do nothing special when writing
            # the HDF5 file (http://www.hdfgroup.org/HDF5/doc/H5.format.html).
            with open(file_name, 'rb') as inp:
                header = inp.read(8)
            return header == HDF5_SIGNATURE
        except Exception:
            return False

    def write(self, h5_file_name=None):
        """
        Write fields to HDF5 file to support restarting the solver
        """
        sim = self.simulation

        if h5_file_name is None:
            h5_file_name = sim.input.get_output_file_path(
                'output/hdf5_file_name', '_savepoint_%08d.h5'
            )
            h5_file_name = h5_file_name % sim.timestep

        # Write dolfin objects using dolfin.HDF5File
        sim.log.info('Creating HDF5 restart file %s' % h5_file_name)
        comm = sim.data['mesh'].mpi_comm()
        with dolfin.HDF5File(comm, h5_file_name, 'w') as h5:
            # Write mesh
            h5.write(sim.data['mesh'], '/mesh')

            # Write mesh facet regions (from mesh generator)
            mfr = sim.data['mesh_facet_regions']
            if mfr is not None:
                h5.write(mfr, '/mesh_facet_regions')

            # Write functions, sorted to ensure deterministic output order
            funcnames = []
            skip = {'coupled'}  # Skip these functions
            funcs = sorted(sim.data.items())
            for name, value in funcs:
                if isinstance(value, dolfin.Function) and name not in skip:
                    assert name not in funcnames
                    h5.write(value, '/%s' % name)
                    funcnames.append(name)

        # Only write metadata on root process
        # Important: no collective operations below this point!
        comm.barrier()
        if self.simulation.rank != 0:
            comm.barrier()
            return h5_file_name

        # Write numpy objects and metadata using h5py.File
        with h5py.File(h5_file_name, 'r+') as hdf:
            # Metadata
            meta = hdf.create_group('ocellaris')
            meta.attrs['time'] = sim.time
            meta.attrs['iteration'] = sim.timestep
            meta.attrs['dt'] = sim.dt
            meta.attrs['restart_file_format'] = 2

            # Functions to save strings
            string_dt = h5py.special_dtype(vlen=str)

            def np_string(root, name, strdata):
                np_data = numpy.array(str(strdata).encode('utf8'), dtype=object)
                root.create_dataset(name, data=np_data, dtype=string_dt)

            def np_stringlist(root, name, strlist):
                np_list = numpy.array(
                    [str(s).encode('utf8') for s in strlist], dtype=object
                )
                root.create_dataset(name, data=np_list, dtype=string_dt)

            # List of names
            repnames = list(sim.reporting.timestep_xy_reports.keys())
            np_stringlist(meta, 'function_names', funcnames)
            np_stringlist(meta, 'report_names', repnames)

            # Save the current input and the full log file
            np_string(meta, 'input_file', sim.input)
            np_string(meta, 'full_log', sim.log.get_full_log())

            # Save reports
            reps = hdf.create_group('reports')
            reps['timesteps'] = numpy.array(sim.reporting.timesteps, dtype=float)
            for rep_name, values in sim.reporting.timestep_xy_reports.items():
                reps[rep_name] = numpy.array(values, dtype=float)

            # Save persistent data dictionaries
            pdd = hdf.create_group('ocellaris_data')
            i = 0
            for name, data in self.persisted_python_data.items():
                # Get stripped down data with only basic data types
                data2 = {}
                for k, v in data.items():
                    if is_basic_datatype(k) and is_basic_datatype(v):
                        data2[k] = v

                if not data2:
                    # no basic data to store, skip this dict
                    continue

                # Convert basic data type data to YAML format
                data3 = yaml.dump(data2)

                pdi = pdd.create_group('data_%02d' % i)
                np_string(pdi, 'name', name)
                np_string(pdi, 'data', data3)
                i += 1

        comm.barrier()
        return h5_file_name

    def read_metadata(self, h5_file_name, function_details=False):
        """
        Read HDF5 restart file metadata
        """
        # Check file format and read metadata
        with h5py.File(h5_file_name, 'r') as hdf:
            if not 'ocellaris' in hdf:
                ocellaris_error(
                    'Error reading restart file',
                    'Restart file %r does not contain Ocellaris meta data'
                    % h5_file_name,
                )

            meta = hdf['ocellaris']
            restart_file_version = meta.attrs['restart_file_format']
            if restart_file_version != 2:
                ocellaris_error(
                    'Error reading restart file',
                    'Restart file version is %d, this version of Ocellaris only '
                    % restart_file_version
                    + 'supports version 2',
                )

            t = float(meta.attrs['time'])
            it = int(meta.attrs['iteration'])
            dt = float(meta.attrs['dt'])
            inpdata = meta['input_file'].value
            funcnames = list(meta['function_names'])

            # Read function signatures
            if function_details:
                signatures = {}
                for fname in funcnames:
                    signatures[fname] = hdf[fname].attrs['signature']

        if function_details:
            return funcnames, signatures
        else:
            return t, it, dt, inpdata, funcnames

    def read(self, h5_file_name, read_input=True, read_results=True):
        """
        Read an HDF5 restart file on the format written by _write_hdf5()
        """
        # Check file format and read metadata
        t, it, dt, inpdata, funcnames = self.read_metadata(h5_file_name)

        sim = self.simulation
        h5 = dolfin.HDF5File(dolfin.MPI.comm_world, h5_file_name, 'r')

        # This flag is used in sim.setup() to to skip mesh creation
        # and may be used by user code etc
        sim.restarted = True

        if read_input:
            # Read the input file
            sim.input.read_yaml(yaml_string=inpdata)
            sim.input.file_name = h5_file_name
            sim.input.set_value('time/tstart', t)
            sim.input.set_value('time/dt', dt)

            # Read mesh data
            mesh = dolfin.Mesh()
            h5.read(mesh, '/mesh', False)
            if h5.has_dataset('/mesh_facet_regions'):
                mesh_facet_regions = dolfin.FacetFunction('size_t', mesh)
                h5.read(mesh_facet_regions, '/mesh_facet_regions')
            else:
                mesh_facet_regions = None
            sim.set_mesh(mesh, mesh_facet_regions)

            # Setup needs to know that the simulation was restarted for
            # XDMF file renaming (among other possible uses)
            sim.timestep = it

        if read_results:
            sim.log.info('Reading fields from restart file %r' % h5_file_name)
            sim.timestep = it

            # Read result field functions
            for name in funcnames:
                sim.log.info('    Function %s' % name)
                h5.read(sim.data[name], '/%s' % name)

            h5.close()  # Close dolfin.HDF5File

            # Read persisted data dictionaries with h5py
            with h5py.File(h5_file_name, 'r') as hdf:
                pdd = hdf.get('ocellaris_data', {})
                for pdi in pdd.values():
                    name = pdi['name'].value
                    data = pdi['data'].value

                    # Parse YAML data and update the persistent data store
                    data2 = yaml.load(data)
                    data3 = self.persisted_python_data.setdefault(name, {})
                    data3.update(data2)

    def read_functions(self, h5_file_name):
        """
        Return a dictionary of Functions read from the HDF5 file
        Mixed or vector function spaces are not supported and
        None will be returned for these
        """
        # Check file format and read metadata
        funcnames, signatures = self.read_metadata(h5_file_name, function_details=True)

        # Read mesh data
        h5 = dolfin.HDF5File(dolfin.MPI.comm_world, h5_file_name, 'r')
        mesh = dolfin.Mesh()
        h5.read(mesh, '/mesh', False)

        def mk_func(name):
            signature = signatures[name].decode('utf8')

            # Parse strings like "FiniteElement('Discontinuous Lagrange', tetrahedron, 2)"
            pattern = r"FiniteElement\('(?P<family>[^']+)', \w+, (?P<degree>\d+)\)"
            m = re.match(pattern, signature)
            if not m:
                return None

            family = m.group('family')
            degree = int(m.group('degree'))

            self.simulation.log.info(
                '        Found function %s of family %r '
                'and degree %d' % (name, family, degree)
            )
            V = dolfin.FunctionSpace(mesh, family, degree)
            return dolfin.Function(V)

        # Read result field functions
        funcs = {}
        for name in funcnames:
            f = mk_func(name)
            if f is not None:
                h5.read(f, '/%s' % name)
            funcs[name] = f
        h5.close()
        return funcs


def is_basic_datatype(value):
    """
    We only save "basic" datatypes like ints, floats, strings and
    lists, dictionaries or sets of these as pickles in restart files

    This function returns True if the datatype is such a basic data
    type
    """
    if isinstance(value, (int, float, str, bytes, bool)):
        return True
    elif isinstance(value, (list, tuple, set)):
        return all(is_basic_datatype(v) for v in value)
    elif isinstance(value, dict):
        return all(
            is_basic_datatype(k) and is_basic_datatype(v) for k, v in value.items()
        )
    return False
