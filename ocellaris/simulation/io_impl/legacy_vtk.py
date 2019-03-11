# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
Write legacy VTK format files (ASCII and BINARY format
See, e.g., http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats

The binary file format writer may be buggy
"""
import numpy
import dolfin
from ocellaris.utils import get_local, ocellaris_error


WRITE_BINARY = False
VTK_QUADRATIC_TRIANGLE = 22
VTK_QUADRATIC_TETRA = 24
UFC2VTK_TET10 = [0, 1, 2, 3, 9, 6, 8, 7, 5, 4]
UFC2VTK_TRI06 = [0, 1, 2, 5, 3, 4]


class LegacyVTKIO:
    def __init__(self, simulation):
        """
        Output to legacy VTK format - Python implementation,
        supports discontinuous quadratic fields
        """
        self.simulation = simulation

    def close(self):
        """
        Save final restart file and close open files
        """
        pass

    def write(self, file_name=None, extra_funcs=(), include_standard=True):
        """
        Write a file that can be used for visualization in Paraview
        """
        sim = self.simulation
        binary_format = sim.input.get_value(
            'output/vtk_binary_format', WRITE_BINARY, 'bool'
        )
        if file_name is None:
            file_name = sim.input.get_output_file_path(
                'output/vtk_file_name', '_%08d.vtk'
            )
            file_name = file_name % sim.timestep

        sim.log.info('    Writing legacy VTK file %s' % file_name)
        with dolfin.Timer('Ocellaris write legacy VTK file'):
            self._write_vtk(file_name, binary_format, extra_funcs, include_standard)

        return file_name

    def _write_vtk(
        self,
        file_name,
        binary_format=WRITE_BINARY,
        extra_funcs=(),
        include_standard=False,
    ):
        """
        Write plot file in legacy VTK format
        """
        sim = self.simulation
        t = numpy.array([sim.time], dtype=float)
        mesh = sim.data['mesh']

        # The functions to output
        funcs = []
        if include_standard:
            for fn in 'u0 u1 u2 rho p p_hydrostatic c'.split():
                func = sim.data.get(fn, None)
                if isinstance(func, dolfin.Function):
                    funcs.append(func)
        funcs.extend(extra_funcs)

        # Get the data from all processes (returns None when not on rank 0)
        results = gather_vtk_info(mesh, funcs)
        if results is not None:
            coords, connectivity, cell_types, func_vals = results
            write_vtk_file(
                file_name, binary_format, coords, connectivity, cell_types, func_vals, t
            )

        mesh.mpi_comm().barrier()


def write_vtk_file(
    file_name, binary_format, coords, connectivity, cell_types, func_vals, t
):
    # VTK writer functions
    def write_ascii(text):
        return out.write(text.encode('ASCII'))

    if binary_format:

        def write_array(data):
            out.write(data.tobytes())
            write_ascii('\n\n')

    else:

        def write_array(data):
            fmt = '%d' if data.dtype == numpy.intc else '%.5E'
            if len(data.shape) == 1:
                write_ascii('\n'.join(fmt % v for v in data))
            else:
                write_ascii('\n'.join(' '.join(fmt % v for v in row) for row in data))
            write_ascii('\n\n')

    # Separate scalars from vectors
    scalars, vectors = {}, {}
    for function_name in list(func_vals.keys()):
        # Vectors are expected to have scalar components named xx0, xx1 and
        # possibly xx2. The xx prefix can be u, up, u_conv, u_star etc
        if not function_name.endswith('0'):
            continue
        prefix = function_name[:-1]
        names = ['%s0' % prefix, '%s1' % prefix, '%s2' % prefix]

        u0 = func_vals.pop(names[0])
        u1 = func_vals.pop(names[1])
        if names[2] in func_vals:
            u2 = func_vals.pop(names[2])
        else:
            u2 = numpy.zeros_like(func_vals['u0'])
        vectors[prefix] = numpy.zeros((len(u0), 3), dtype=float)
        vectors[prefix][:, 0] = u0
        vectors[prefix][:, 1] = u1
        vectors[prefix][:, 2] = u2
    scalars.update(func_vals)

    # Write the VTK file
    with open(file_name, 'wb') as out:
        # Write geometry
        write_ascii('# vtk DataFile Version 3.0\n')
        write_ascii('Ocellaris simulation output\n')
        if binary_format:
            write_ascii('BINARY\n')
        else:
            write_ascii('ASCII\n')

        Nverts = coords.shape[0]
        Ncells = connectivity.shape[0]

        write_ascii('DATASET UNSTRUCTURED_GRID\n')
        write_ascii('FIELD FieldData 1\n')
        write_ascii('TIME 1 1 double\n')  # Supported in VisIt, not Paraview
        write_array(t)

        write_ascii('POINTS %d float\n' % Nverts)
        write_array(coords)

        write_ascii('CELLS %d %d\n' % (Ncells, Ncells * 11))
        write_array(connectivity)

        write_ascii('CELL_TYPES %d\n' % Ncells)
        write_array(cell_types)

        write_ascii('POINT_DATA %d\n' % Nverts)

        # Write scalar functions
        for function_name in sorted(scalars):
            write_ascii('SCALARS %s float 1\n' % function_name)
            write_ascii('LOOKUP_TABLE default\n')
            write_array(scalars[function_name])

        # Write vector functions
        for function_name in sorted(vectors):
            write_ascii('VECTORS %s float\n' % function_name)
            write_array(vectors[function_name])


def gather_vtk_info(mesh, funcs):
    """
    Gather the necessary information to write a legacy VTK output file
    on the root process
    """
    # The code below assumes that the first function is DG2 (u0)
    assert funcs[0].function_space().ufl_element().family() in (
        'Discontinuous Lagrange',
        'Lagrange',
    )
    assert funcs[0].function_space().ufl_element().degree() == 2

    # The code is currently 3D only
    gdim = mesh.geometry().dim()
    dofs_x = funcs[0].function_space().tabulate_dof_coordinates().reshape((-1, gdim))
    assert gdim in (2, 3), 'VTK output currently only supported for 2D and 3D meshes'

    # Collect information about the functions
    func_names = []
    all_vals = []
    dofmaps = []
    for u in funcs:
        func_names.append(u.name())
        all_vals.append(get_local(u))
        dofmaps.append(u.function_space().dofmap())

    # Collect local data
    local_res = _collect_3D(mesh, gdim, dofs_x, func_names, all_vals, dofmaps)

    # MPI communication to get all data on root process
    comm = mesh.mpi_comm()
    all_results = comm.gather(local_res)
    if all_results is None:
        return None

    # Assemble information from all processes
    coords, connectivity, cell_types = [], [], []
    func_vals = {n: [] for n in func_names}
    for coords_i, connectivity_i, cell_types_i, func_vals_i in all_results:
        K = len(coords)
        coords.extend(coords_i)
        for conn in connectivity_i:
            connectivity.append(conn[:1] + [c + K for c in conn[1:]])
        cell_types.extend(cell_types_i)
        for n in func_names:
            func_vals[n].extend(func_vals_i[n])

    # Convert to numpy arrays
    coords = numpy.array(coords, dtype=numpy.float32)
    connectivity = numpy.array(connectivity, dtype=numpy.intc)
    cell_types = numpy.array(cell_types, dtype=numpy.intc)
    for n in func_names:
        func_vals[n] = numpy.array(func_vals[n], dtype=numpy.float32)

    # Check that the data is appropriate for output
    Nverts = len(coords)
    Ncells = len(connectivity)
    assert coords.shape == (Nverts, 3)
    assert connectivity.shape == (Ncells, 11 if gdim == 3 else 7)
    assert cell_types.shape == (Ncells,)
    for n in func_names:
        assert func_vals[n].shape == (Nverts,)

    return coords, connectivity, cell_types, func_vals


def _collect_3D(mesh, gdim, dofs_x, func_names, all_vals, dofmaps):
    """
    Collect geometry and function values on local process
    """
    if gdim == 3:
        # Assume max degree is 2
        VTK_TYPE = VTK_QUADRATIC_TETRA
        UFC2VTK = UFC2VTK_TET10

    elif gdim == 2:
        # Assume max degree is 2
        VTK_TYPE = VTK_QUADRATIC_TRIANGLE
        UFC2VTK = UFC2VTK_TRI06

        # Add z-coordinate
        dx = numpy.zeros((len(dofs_x), 3), float)
        dx[:, :2] = dofs_x
        dofs_x = dx

    coords, connectivity, cell_types = [], [], []
    func_vals = {n: [] for n in func_names}
    M = len(UFC2VTK)
    for cell in dolfin.cells(mesh, 'regular'):
        cidx = cell.index()

        # Get the cell geometry and connectivity
        # (assumes that the max degree is 2)
        cell_types.append(VTK_TYPE)
        connectivity.append([M])
        dofs = dofmaps[0].cell_dofs(cidx)
        Lc = len(coords)
        for k in range(M):
            d = dofs[k]
            coords.append(tuple(dofs_x[d]))
            connectivity[-1].append(Lc + UFC2VTK[k])

        for name, vals, dm in zip(func_names, all_vals, dofmaps):
            dofs = dm.cell_dofs(cidx)
            m = len(dofs)

            if gdim == 3 and m == 10:
                cvals = vals[dofs]

            elif gdim == 3 and m == 4:
                cvals = list(vals[dofs])
                cvals.append((cvals[2] + cvals[3]) / 2)
                cvals.append((cvals[1] + cvals[3]) / 2)
                cvals.append((cvals[1] + cvals[2]) / 2)
                cvals.append((cvals[0] + cvals[3]) / 2)
                cvals.append((cvals[0] + cvals[2]) / 2)
                cvals.append((cvals[0] + cvals[1]) / 2)

            elif gdim == 3 and m == 1:
                cvals = [vals[dofs[0]]] * M

            elif gdim == 2 and m == 6:
                cvals = vals[dofs]

            elif gdim == 2 and m == 3:
                cvals = list(vals[dofs])
                cvals.append((cvals[1] + cvals[2]) / 2)
                cvals.append((cvals[0] + cvals[2]) / 2)
                cvals.append((cvals[0] + cvals[1]) / 2)

            elif gdim == 2 and m == 1:
                cvals = [vals[dofs[0]]] * M

            else:
                ocellaris_error(
                    'VTK write error',
                    'Unsupported geometry dimension / ' 'element type for %s' % name,
                )
            func_vals[name].extend(cvals)

    return coords, connectivity, cell_types, func_vals
