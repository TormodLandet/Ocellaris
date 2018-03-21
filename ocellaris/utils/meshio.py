import os, uuid
import dolfin
from ocellaris.utils import ocellaris_error


def load_meshio_mesh(mesh, file_name, meshio_file_type=None):
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
    
    # Open mesh for editing
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


def build_distributed_mesh(mesh):
    """
    Work around missing dolfin pybind11 wrapped method
    FIXME: get this into dolfin / dolfin-x
    """
    if build_distributed_mesh.func is None:
        cpp_code = """
        #include<pybind11/pybind11.h>
        #include<dolfin/mesh/MeshPartitioning.h>
        #include<dolfin/mesh/Mesh.h>

        void my_distributer(dolfin::Mesh mesh)
        {
            return dolfin::MeshPartitioning::build_distributed_mesh(mesh);
        }

        namespace py = pybind11;

        PYBIND11_MODULE(SIGNATURE, m) {
           m.def("build_distributed_mesh", &my_distributer, py::arg("mesh"));
        }
        """
        build_distributed_mesh.func = dolfin.compile_cpp_code(cpp_code).build_distributed_mesh
    return build_distributed_mesh.func(mesh)
build_distributed_mesh.func = None

