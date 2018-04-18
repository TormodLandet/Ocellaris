import os
import dolfin
import numpy
from collections import OrderedDict
from ocellaris.utils import init_mesh_geometry, timeit, ocellaris_error
from . import Probe, register_probe


WRITE_INTERVAL = 1


@register_probe('PlaneProbe')
class PlaneProbe(Probe):
    description = 'Produce a 2D slice function from a 3D mesh'
    
    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        assert self.simulation.ndim == 3, 'PlaneProbe only implemented in 3D'

        # Read input
        inp = probe_input
        self.name = inp.get_value('name', required_type='string')
        self.plane_point = inp.get_value('plane_point', required_type='list(float)',
                                         required_length=3)
        self.plane_normal = inp.get_value('plane_normal', required_type='list(float)',
                                          required_length=3)
        self.custom_hook_point = inp.get_value('custom_hook', None, required_type='string')
        self.write_interval = inp.get_value('write_interval', WRITE_INTERVAL, 'int')
        
        # Get the names of the function(s) to be sliced
        fn = inp.get_value('field', required_type='any')
        if isinstance(fn, str):
            self.field_names = [fn]
        else:
            self.field_names = inp.validate_and_convert('field', fn, 'list(string)')
        
        # Get the functions and verify the function spaces
        self.funcs_3d = []
        for fn in self.field_names:
            func_3d = simulation.data[fn]
            V = func_3d.function_space()
            fam = V.ufl_element().family()
            deg = V.ufl_element().degree()
            if not self.funcs_3d:
                self.family = fam
                self.degree = deg
            elif fam != self.family or deg != self.degree:
                ocellaris_error('Mismatching function spaces in PlainProbe %s' % self.name,
                                'All functions must have the same function space. ' +
                                '%s is %r but %r was expected' % (fn, (fam, deg),
                                                                  (self.family,
                                                                   self.degree)))
            self.funcs_3d.append(func_3d)
        
        # Create the slice
        self.slice = FunctionSlice(self.plane_point, self.plane_normal, V)
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        self.file_name = '%s_slice_%s.xdmf' % (prefix, self.name)
        
        if simulation.rank == 0:
            V_2d = self.slice.slice_function_space
            mesh_2d = V_2d.mesh()
            simulation.log.info('        Created 2D mesh with %r cells' % mesh_2d.num_cells())
            
            # Remove any previous XDMF files
            file_name = self.file_name
            file_name2 = os.path.splitext(file_name)[0] + '.h5'
            if os.path.isfile(file_name):
                simulation.log.info('        Removing existing XDMF file %s' % file_name)
                os.remove(file_name)
            if os.path.isfile(file_name2):
                simulation.log.info('        Removing existing XDMF file %s' % file_name2)
                os.remove(file_name2)
            
            simulation.log.info('        Creating XDMF file %s' % self.file_name)
            self.xdmf_file = dolfin.XDMFFile(dolfin.MPI.comm_self, file_name)
            self.xdmf_file.parameters['flush_output'] = True
            self.xdmf_file.parameters['rewrite_function_mesh'] = False
            self.xdmf_file.parameters['functions_share_mesh'] = True
        
            # Create storage for 2D functions
            self.funcs_2d = []
            for fn in self.field_names:
                func_2d = dolfin.Function(V_2d)
                func_2d.rename(fn, fn)
                self.funcs_2d.append(func_2d)
        else:
            self.funcs_2d = [None] * len(self.field_names)
        
        if self.custom_hook_point is not None:
            simulation.hooks.add_custom_hook(self.custom_hook_point, self.run, 'Probe "%s"' % self.name)
        else:
            self.end_of_timestep = self.run
    
    def run(self):
        """
        Find and output the plane probe
        """
        it = self.simulation.timestep
        
        # Should we update the file?
        if not (it == 1 or it % self.write_interval == 0):
            return
        
        for func_3d, func_2d in zip(self.funcs_3d, self.funcs_2d):
            self.slice.get_slice(func_3d, func_2d)
            if self.simulation.rank == 0:
                self.xdmf_file.write(func_2d, self.simulation.time)


class FunctionSlice:
    def __init__(self, pt, n, V3d):
        """
        Take the definition of a plane and a 3D function space
        Construct a 2D mesh on rank 0 (only) and the necessary
        data structures to extract function values at the 2D
        mesh in an efficient way
        
        * pt: a point in the plane
        * n: a normal vector to the plane. Does not need to be a unit normal
        * V3d: the 3D function space to be intersected by the plane 
        """
        gdim = V3d.mesh().geometry().dim()
        assert gdim == 3, 'Function slice only supported in 3D'
        
        # 3D function space data
        comm = V3d.mesh().mpi_comm()
        elem_3d = V3d.ufl_element()
        family = elem_3d.family()
        degree = elem_3d.degree()
        
        # Create the 2D mesh
        # The 2D mesh uses MPI_COMM_SELF and lives only on the root process
        mesh_2d, cell_origins = make_cut_plane_mesh(pt, n, V3d.mesh())
        
        # Precompute data on root process
        if comm.rank == 0:
            # Make the 2D function space
            
            V2d = dolfin.FunctionSpace(mesh_2d, family, degree)
            self.slice_function_space = V2d
            
            # Get 2D dof coordinates and dofmap
            dof_pos_2d = V2d.tabulate_dof_coordinates().reshape((-1, gdim))
            dofmap_2d = V2d.dofmap()
            
            # Link 2D dofs and 3D ranks
            links_for_rank = [[] for _ in range(comm.size)]
            for cell in dolfin.cells(mesh_2d):
                cid = cell.index()
                
                # Assume no cell renumbering
                orig_rank, orig_cell_index = cell_origins[cid]
                
                for dof in dofmap_2d.cell_dofs(cid):
                    links_for_rank[orig_rank].append((dof, orig_cell_index))
            
            # Distribute data to all ranks
            distribute_this = []
            for rank in range(comm.size):
                positions = [dof_pos_2d[i] for i, _ in links_for_rank[rank]]
                orig_cells = [ocid for _, ocid in links_for_rank[rank]]
                positions = numpy.array(positions, float)
                orig_cells = numpy.array(orig_cells, numpy.intc)
                distribute_this.append((positions, orig_cells))
            
            # Store which 2D dof belongs on which rank
            self._dofs_for_rank = []
            for rank in range(comm.size):
                dfr = [dof for dof, _ in links_for_rank[rank]]
                self._dofs_for_rank.append(numpy.array(dfr, int))
        else:
            distribute_this = None
        
        # Get positions along with the index of the 3D cell for all points that
        # need to be evaluated in order to build the 2D function
        # Each rank gets positions corresponding to cells located on that rank
        positions, cell_index_3d = comm.scatter(distribute_this)
        
        # Establish efficient ways to get the 2D data from the 3D function
        cell_dofs = [V3d.dofmap().cell_dofs(i) for i in cell_index_3d]
        self._cell_dofs = numpy.array(cell_dofs, int)
        self._factors = numpy.zeros(self._cell_dofs.shape, float)
        self._local_data = numpy.zeros(len(cell_dofs), float)
        evaluate_basis_functions(V3d, positions, cell_index_3d, self._factors)
    
    @timeit.named('FunctionSlice.get_slice')
    def get_slice(self, func_3d, func_2d=None):
        """
        Return the function on the 2D slice of the 3D mesh 
        """
        comm = func_3d.function_space().mesh().mpi_comm()
        
        # Get local values from all processes
        arr_3d = func_3d.vector().get_local()
        local_data = self._local_data
        N = local_data.size
        facs = self._factors
        cd = self._cell_dofs
        for i in range(N):
            local_data[i] = arr_3d[cd[i]].dot(facs[i]) 
        all_data = comm.gather(local_data)
        
        if comm.rank == 0:
            if func_2d is None:
                func_2d = dolfin.Function(self.slice_function_space)
            arr_2d = func_2d.vector().get_local()
            for data, dofs in zip(all_data, self._dofs_for_rank):
                arr_2d[dofs] = data
            func_2d.vector().set_local(arr_2d)
            func_2d.vector().apply('insert')
            
            return func_2d 


def make_cut_plane_mesh(pt, n, mesh3d):
    """
    Returns a 2D mesh of the intersection of a 3D mesh and a plane
    
    * pt: a point in the plane
    * n: a normal vector to the plane. Does not need to be a unit normal
    * mesh3d: the 3D mesh to be intersected by the plane
    
    This function assumes that the 3D mesh consists solely of tetrahedra which
    gives a 2D mesh of triangles 
    """
    # Get results on this rank and send to root process
    rank_results = get_points_in_plane(pt, n, mesh3d)
    comm = mesh3d.mpi_comm()
    all_results = comm.gather((comm.rank, rank_results))
    
    # No collective operatione below this point!
    if comm.rank != 0:
        return None, None
    
    point_ids = {}
    points = []
    connectivity = []
    cell_origins = []
    for rank, res in all_results:
        for cell_id, cell_coords in res.items():
            if len(cell_coords) == 4:
                subcells = (cell_coords[1:], cell_coords[:-1])
            else:
                subcells = (cell_coords,)
            
            for cell_coords in subcells:
                cell_points = []
                for coords in cell_coords: 
                    if coords not in point_ids:
                        point_ids[coords] = len(point_ids)
                        points.append(coords)
                    cell_points.append(point_ids[coords])
                connectivity.append(cell_points)
                cell_origins.append((rank, cell_id))
    
    # Create the mesh
    points = numpy.array(points, float)
    connectivity = numpy.array(connectivity, int)
    tdim, gdim = 2, 3
    mesh2d = dolfin.Mesh(dolfin.MPI.comm_self)
    init_mesh_geometry(mesh2d, points, connectivity, tdim, gdim)
    
    return mesh2d, cell_origins


def get_points_in_plane(pt, n, mesh, eps=1e-8):
    """
    Returns a dictionary of cell_id -> list of three/four coordinate (x, y, z)
    tuples for any cell that has three points on a mesh edge that crosses the
    given plane. For edges coincident with the plane only the two end points
    are returned
    
    * pt, n: needed input to get_plane_coefficients()
    * mesh: only the local part of the mesh is explored, no MPI communication
    * eps: distance to be considered in the plane or not
    """
    assert mesh.geometry().dim() == 3, 'Can only find planes in 3D meshes'
    conn_C_V = mesh.topology()(3, 0)
    coords = mesh.coordinates()
    
    plane_coefficients = get_plane_coefficients(pt, n)
    all_sides = get_point_sides(plane_coefficients, coords)
    sign = lambda x: -1 if x < 0 else 1
    n = plane_coefficients[:3] # normal to the plane
    
    results = OrderedDict()
    for cell in dolfin.cells(mesh, 'regular'):
        cid = cell.index()
        verts = conn_C_V(cid)
        sides = [all_sides[v] for v in verts]
        cell_points = []
        
        # Check for points on plane
        on_plane = [abs(s) < eps for s in sides]
        for i, v in enumerate(verts):
            if on_plane[i]:
                pos = coords[v]
                cell_points.append((pos[0], pos[1], pos[2]))
        
        for v0, v1 in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
            if on_plane[v0] or on_plane[v1]:
                continue
            elif sign(sides[v0]) == sign(sides[v1]):
                continue
            
            # Points v0 and v1 are on different sides of the plane
            # Get the mesh-local vertex numbers in ascending order 
            vid0, vid1 = verts[v0], verts[v1]
            if vid1 < vid0:
                vid1, vid0 = vid0, vid1
            
            # Get the crossing point
            c0 = coords[vid0]
            c1 = coords[vid1]
            u = c1 - c0
            f = numpy.dot(pt - c0, n) / numpy.dot(u, n)
            pos = c0 + f * u
            cell_points.append((pos[0], pos[1], pos[2]))
        
        # Only add this cell if in contributes with some area to the plane
        if len(cell_points) > 2:
            results[cid] = cell_points
    
    return results


def get_plane_normal(p1, p2, p3):
    """
    Given a plane passing through the points p1, p2 and p3
    find the normal vector to this plane
    """
    p1 = numpy.array(p1, float)
    p2 = numpy.array(p2, float)
    p3 = numpy.array(p3, float)
    
    u = p1 - p2
    v = p1 - p3
    w = numpy.cross(u, v)
    
    assert u.dot(u) > 1e-6, 'Points p1 and p2 coincide!'
    assert v.dot(v) > 1e-6, 'Points p1 and p3 coincide!'
    assert w.dot(w) > 1e-6, 'Points p1, p2 and p3 fall on the same line!'
    
    return w / (w.dot(w))**0.5


def get_plane_coefficients(pt, n):
    """
    Get the equation for the plane containing point pt when n is the normal
    vector. Returns A, B, C, D from the function
    
        f(x, y, z) = A x + B y + C z + D
    
    which can be evaluated at any point and is zero on the plane, -1 on one
    side and +1 on the other side of the plane 
    """
    
    return n[0], n[1], n[2], -numpy.dot(pt, n)


def get_point_sides(plane_coeffs, points):
    """
    Given the four plane coefficients (as computed by get_plane_coefficients)
    and a numpy array of shape (3, N) with N points return a boolean array
    of N values with numbers that are positive or negative depending which side
    of the plane the corresponding point is located. Values (close to) zero are
    coincident with the plane
    """
    return numpy.dot(points, plane_coeffs[:3]) + plane_coeffs[3]


def evaluate_basis_functions(V, positions, cell_indices, factors):
    """
    Current FEniCS pybind11 bindings lack wrappers for these functions,
    so a small C++ snippet is used to generate the necessary dof factors
    to evaluate a function at the given points
    """
    if evaluate_basis_functions.func is None:
        cpp_code = """
        #include <vector>
        #include <pybind11/pybind11.h>
        #include <pybind11/eigen.h>
        #include <pybind11/numpy.h>
        #include <dolfin/fem/FiniteElement.h>
        #include <dolfin/function/FunctionSpace.h>
        #include <dolfin/mesh/Mesh.h>
        #include <dolfin/mesh/Cell.h>
        #include <Eigen/Core>
        
        using IntVecIn = Eigen::Ref<const Eigen::VectorXi>;
        using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;        
        namespace py = pybind11;
        
        void dof_factors(const dolfin::FunctionSpace &V, const RowMatrixXd &positions,
                         const IntVecIn &cell_indices, Eigen::Ref<RowMatrixXd> out)
        {
            const int N = out.rows();
            if (N == 0)
                return;
            
            const auto &element = V.element();
            const auto &mesh = V.mesh();
            const auto &ufc_element = element->ufc_element();
            std::vector<double> coordinate_dofs;
            
            const std::size_t size = ufc_element->value_size();
            const std::size_t space_dimension = ufc_element->space_dimension();
            if (size * space_dimension != out.cols())
                throw std::length_error("ERROR: out.cols() != ufc element size * ufc element space_dimension");
            
            for (int i = 0; i < N; i++)
            {
                int cell_index = cell_indices(i);
                dolfin::Cell cell(*mesh, cell_index);
                cell.get_coordinate_dofs(coordinate_dofs);
                element->evaluate_basis_all(out.row(i).data(),
                                            positions.row(i).data(),
                                            coordinate_dofs.data(),
                                            cell.orientation());
            }
        }
        
        PYBIND11_MODULE(SIGNATURE, m) {
           m.def("dof_factors", &dof_factors, py::arg("V"),
                 py::arg("positions"), py::arg("cell_indices"),
                 py::arg("out").noconvert());
        }
        """
        evaluate_basis_functions.func = dolfin.compile_cpp_code(cpp_code).dof_factors
    
    assert len(cell_indices) == len(positions) == len(factors)
    evaluate_basis_functions.func(V._cpp_object, positions, cell_indices, factors)
evaluate_basis_functions.func = None
