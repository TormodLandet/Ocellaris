import os, uuid
import dolfin
from ocellaris.utils import ocellaris_error


def load_meshio_mesh(mpi_comm, file_name, meshio_file_type=None):
    """
    Use Nico Schlömer's meshio library to read a mesh
    https://github.com/nschloe/meshio
    
    The meshio_file_type can be None (auto detect based on file extension) or
    one of the following (taken from meshio source code in March 2018):
    
    * ansys
    * exodus
    * gmsh-ascii
    * gmsh-binary
    * dolfin-xml
    * med
    * medit
    * permas
    * moab
    * off
    * stl-ascii
    * stl-binary
    * vtk-ascii
    * vtk-binary
    * vtu-ascii
    * vtu-binary
    * xdmf
    
    """
    # Construct mesh on rank 0
    if mpi_comm.rank == 0:
        # Create a h5 file containing the mesh by use of meshio
        comm_self = dolfin.MPI.comm_self 
        mesh, _ = _make_mesh(comm_self, file_name, meshio_file_type)
        
        # Find an available file name
        while True:
            h5_file_name = '%s.%s.h5' % (file_name, uuid.uuid4())
            if not os.path.exists(h5_file_name):
                break
        
        # Save to h5 format
        with dolfin.HDF5File(comm_self, h5_file_name, "w") as h5:
            h5.write(mesh, '/mesh')
    
    # Get the file name to all processes
    h5_file_name = mpi_comm.bcast(h5_file_name, root=0)
    
    # Load the mesh
    mesh = dolfin.Mesh(mpi_comm)
    with dolfin.HDF5File(mpi_comm, h5_file_name, "r") as h5:
        h5.read(mesh, '/mesh', False)
    
    mpi_comm.barrier()
    if mpi_comm.rank == 0:
        os.remove(h5_file_name)
    
    # Return mesh and facet markers
    return mesh, None


def _make_mesh(mpi_comm, file_name, meshio_file_type=None):
    try:
        import meshio
    except ImportError:
        ocellaris_error('Please install meshio',
                        'The meshio library is missing. Run something like '
                        '"python3 -m pip install --user meshio" to install it.')
    
    points, cells, _point_data, _cell_data, _field_data = \
        meshio.read(file_name, meshio_file_type)
    
    if 'triangle' in cells:
        dim = 2
        cell_type = 'triangle'
        connectivity = cells['triangle']
    elif 'tetra' in cells:
        dim = 3
        cell_type = 'tetrahedron'
        connectivity = cells['tetra']
    else:
        ocellaris_error('Mesh loaded by meshio contains unsupported cells',
                        'Expected to find tetrahedra or triangles, found %r'
                        % (tuple(cells.keys()), ))
    assert 'triangle' in cells or 'tetra' in cells
    
    # Create the mesh and open for editing
    mesh = dolfin.Mesh(mpi_comm)
    editor = dolfin.MeshEditor()
    editor.open(mesh, cell_type, dim, dim)
    
    # Add vertices
    editor.init_vertices(points.shape[0])
    for i, p in enumerate(points):
        editor.add_vertex(i, p)
    
    # Add cells
    editor.init_cells(connectivity.shape[0])
    for i, c in enumerate(connectivity):
        editor.add_cell(i, c)
    
    editor.close()
    return mesh, None
