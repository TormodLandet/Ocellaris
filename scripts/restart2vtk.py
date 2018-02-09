"""
Take one function from an Ocellaris restart h5 file and
export it as a true DG field to a *.vtk file

Currently only implemented for scalar DG2 fields, should
be easy to extend to other element types. The binary file
implementation may be buggy

See, e.g., http://www.earthmodels.org/software/vtk-and-paraview/vtk-file-formats
"""
import sys
import re
import time
from contextlib import contextmanager
import numpy
import h5py
import dolfin


VTK_QUADRATIC_TETRA = 24
UFC2VTK_TET10 = [0, 1, 2, 3, 9, 6, 8, 7, 5, 4]

@contextmanager
def timer(name):
    t1 = time.time()
    print(name, '...')
    yield
    d = time.time() - t1
    print(name, 'DONE in %.2fs' % d)


def read_h5_function(h5_file_name, func_name):
    with timer('  Read function signature'):
        with h5py.File(h5_file_name, 'r') as hdf:
            signature = hdf[func_name].attrs['signature'].decode('utf8')
        
        # Parse strings like "FiniteElement('Discontinuous Lagrange', tetrahedron, 2)"
        pattern = r"FiniteElement\('(?P<family>[^']+)', \w+, (?P<degree>\d+)\)"
        m = re.match(pattern, signature)
        if not m:
            return None
        family = m.group('family')
        degree = int(m.group('degree'))
        print('    Got', family, degree)
    
    # Read the mesh and the function
    with dolfin.HDF5File(dolfin.MPI.comm_world, h5_file_name, 'r') as h5:
        with timer('  Read mesh'):
            mesh = dolfin.Mesh()
            h5.read(mesh, '/mesh', False)
        
        with timer('  Make V'):
            V = dolfin.FunctionSpace(mesh, family, degree)
        
        with timer('  Make u'):
            u = dolfin.Function(V)
        
        with timer('  Read u'):
            h5.read(u, '/%s' % func_name)
    return u


def get_geom_and_dof_data(u):
    V = u.function_space()
    dm = V.dofmap()
    mesh = V.mesh()
    gdim = mesh.geometry().dim()
    assert gdim == 3
    dofs_x = V.tabulate_dof_coordinates().reshape((-1, gdim))
    vals = u.vector().get_local()
    
    coords = []
    connectivity = []
    dof_vals = []
    cell_types = []
    for i, cell in enumerate(dolfin.cells(mesh)):
        if i % 1000 == 0 and i > 0:
            print(i)
        
        dofs = dm.cell_dofs(cell.index())
        M = len(dofs)
        cell_types.append(VTK_QUADRATIC_TETRA)
        j = len(coords)
        
        connectivity.append([M])
        for k in range(M):
            d = dofs[k]
            coords.append(tuple(dofs_x[d]))
            connectivity[-1].append(j + UFC2VTK_TET10[k])
            dof_vals.append(vals[d])
    
    return coords, cell_types, connectivity, dof_vals


def restart_file_to_vtk(h5_file_name, function_name, vtk_file_name, binary=False):
    with timer('Reading h5 file %s' % h5_file_name):
        u = read_h5_function(h5_file_name, function_name)
    
    with timer('Extracting geometry and the %s field' % function_name):
        coords, cell_types, connectivity, dof_vals = get_geom_and_dof_data(u)
        
        Nverts = len(coords)
        Ncells = len(connectivity)
        coords = numpy.array(coords, dtype=numpy.float32)
        connectivity = numpy.array(connectivity, dtype=numpy.intc)
        cell_types = numpy.array(cell_types, dtype=numpy.intc)
        dof_vals = numpy.array(dof_vals, dtype=numpy.float32)
        
        for arrname in 'coords connectivity cell_types dof_vals'.split():
            a = locals()[arrname]
            print('    %s has shape %r and dtype %r' % (arrname, a.shape, a.dtype)) 
        
        assert coords.shape == (Nverts, 3)
        assert connectivity.shape == (Ncells, 11)
        assert cell_types.shape == (Ncells,)
        assert dof_vals.shape == (Nverts,)
    
    with timer('Writing VTK file %s' % vtk_file_name), open(vtk_file_name, 'wb') as out:
        write_ascii = lambda text: out.write(text.encode('ASCII'))
        
        write_ascii('# vtk DataFile Version 3.0\n')
        write_ascii('Ocellaris field %s from %s\n' % (function_name, h5_file_name[-100:]))
        
        if binary:
            write_ascii('BINARY\n')
            def write_array(data):
                out.write(data.tobytes())
                write_ascii('\n\n')
        else:
            write_ascii('ASCII\n')
            def write_array(data):
                fmt = '%d' if data.dtype == numpy.intc else '%.5E'
                if len(data.shape) == 1:
                    write_ascii(' '.join(fmt % v for v in data))
                else:
                    write_ascii('\n'.join(' '.join(fmt % v for v in row) for row in data))
                write_ascii('\n\n')
        
        write_ascii('DATASET UNSTRUCTURED_GRID\n')
        write_ascii('POINTS %d float\n' % Nverts)
        write_array(coords)
        
        write_ascii('CELLS %d %d\n' % (Ncells, Ncells*11))
        write_array(connectivity)
        
        write_ascii('CELL_TYPES %d\n' % Ncells)
        write_array(cell_types)
        
        write_ascii('POINT_DATA %d\n' % Nverts)
        
        write_ascii('SCALARS %s float 1\n' % function_name)
        write_ascii('LOOKUP_TABLE default\n')
        write_array(dof_vals)
    
    print('DONE')


if __name__ == '__main__':
    h5_file_name = sys.argv[1]
    function_name = sys.argv[2]
    vtk_file_name = h5_file_name + '.vtk'
    restart_file_to_vtk(h5_file_name, function_name, vtk_file_name)
