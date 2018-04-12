import dolfin
from ocellaris.utils import ocellaris_error


def load_meshio_mesh(mesh, file_name, meshio_file_type=None, sort_order=None):
    """
    Use Nico Schlömer's meshio library to read a mesh, see 
    https://github.com/nschloe/meshio for more information on meshio.
    
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
    
    The sort_order order keyword, if specified, must be an iterable of
    integers. Giving sort_order=(2, 0, 1) will sort elements first by
    z-coordinate, then by x and last by y. Vertices are not sorted.
    """
    try:
        import meshio
    except ImportError:
        ocellaris_error('Please install meshio',
                        'The meshio library is missing. Run something like '
                        '"python3 -m pip install --user meshio" to install it.')
    
    points, cells, _point_data, _cell_data, _field_data = \
        meshio.read(file_name, meshio_file_type)
    
    if 'triangle' in cells:
        assert len(cells) == 1, 'Mixed element meshes is not supported'
        dim = 2
        connectivity = cells['triangle']
    elif 'tetra' in cells:
        assert len(cells) == 1, 'Mixed element meshes is not supported'
        dim = 3
        connectivity = cells['tetra']
    else:
        ocellaris_error('Mesh loaded by meshio contains unsupported cells',
                        'Expected to find tetrahedra or triangles, found %r'
                        % (tuple(cells.keys()), ))
    
    # Order elements by location of the first vertex
    if sort_order is not None:
        connectivity = [[int(c) for c in conn] for conn in connectivity]
        connectivity.sort(key=lambda el: tuple(points[el[0]][i] for i in sort_order)) 
    
    # Add the vertices and cells
    init_mesh_geometry(mesh, points, connectivity, dim, dim)


def init_mesh_geometry(mesh, points, connectivity, tdim, gdim):
    """
    Create a dolfin mesh from a list of points and connectivity  for each cell
    (as returned by meshio). The geometric dimmension dim should be 2 or 3
    """
    cell_type = {2: 'triangle', 3: 'tetrahedron'}[tdim]
    
    # Get global mesh sizes
    Nvert = points.shape[0]
    Ncell = len(connectivity)
    
    # Open mesh for editing
    editor = dolfin.MeshEditor()
    editor.open(mesh, cell_type, tdim, gdim)
    
    # Add vertices
    editor.init_vertices_global(Nvert, Nvert)
    for i, pnt in enumerate(points):
        editor.add_vertex(i, [float(xi) for xi in pnt])
    
    # Add cells
    editor.init_cells_global(Ncell, Ncell)
    for i, conn in enumerate(connectivity):
        editor.add_cell(i, conn)
    
    editor.close()


def build_distributed_mesh(mesh):
    """
    Work around some missing dolfin pybind11 wrapped methods
    Performs the same actions as what happens when reading
    a mesh from an XML file, except for the very first step
    (constructing the whole global mesh on the root process)
    which is assumed to have happened allready 
    """
    if build_distributed_mesh.func is None:
        cpp_code = """
        #include <pybind11/pybind11.h>
        #include <dolfin/common/MPI.h>
        #include <dolfin/mesh/MeshPartitioning.h>
        #include <dolfin/mesh/LocalMeshData.h>
        #include <dolfin/mesh/Mesh.h>
        #include <dolfin/parameter/GlobalParameters.h>
        
        /**
            Code taken from XMLFile::read(Mesh& input_mesh)
            in dolfin/io/XMLFile.cpp on 2018-03-21
        */
        void distribute_mesh(dolfin::Mesh &mesh)
        {
            // Distribute local data
            mesh.domains().clear();
            dolfin::LocalMeshData local_mesh_data(mesh);
            
            // Mesh domain data not present
            if (dolfin::MPI::rank(mesh.mpi_comm()) == 0)
                local_mesh_data.domain_data.clear();
            
            // Partition and build mesh
            const std::string ghost_mode = dolfin::parameters["ghost_mode"];
            dolfin::MeshPartitioning::build_distributed_mesh(mesh, local_mesh_data, ghost_mode);
        }
        
        namespace py = pybind11;
        
        PYBIND11_MODULE(SIGNATURE, m) {
           m.def("distribute_mesh", &distribute_mesh, py::arg("mesh"));
        }
        """
        build_distributed_mesh.func = dolfin.compile_cpp_code(cpp_code).distribute_mesh
    return build_distributed_mesh.func(mesh)
build_distributed_mesh.func = None

