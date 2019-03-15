# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
import numpy
from collections import OrderedDict
from ocellaris.utils import init_mesh_geometry, timeit, ocellaris_error
from ocellaris.simulation.io_impl.xdmf import get_xdmf_file_name
from . import Probe, register_probe


WRITE_INTERVAL = 1


@register_probe('PlaneProbe')
class PlaneProbe(Probe):
    description = 'Produce a 2D slice function from a 3D mesh'

    def __init__(self, simulation, probe_input):
        self.simulation = simulation
        self.family = None
        self.degree = None

        assert self.simulation.ndim == 3, 'PlaneProbe only implemented in 3D'
        assert not self.simulation.mesh_morpher.active, 'PlaneProbe does not support ALE yet'

        # Read input
        inp = probe_input
        self.name = inp.get_value('name', required_type='string')
        self.plane_point = inp.get_value(
            'plane_point', required_type='list(float)', required_length=3
        )
        self.plane_normal = inp.get_value(
            'plane_normal', required_type='list(float)', required_length=3
        )
        xlim = inp.get_value('xlim', None, 'list(float)', required_length=2)
        ylim = inp.get_value('ylim', None, 'list(float)', required_length=2)
        zlim = inp.get_value('zlim', None, 'list(float)', required_length=2)

        # Get the names of the function(s) to be sliced
        fn = inp.get_value('field', required_type='any')
        if isinstance(fn, str):
            self.field_names = [fn]
        else:
            self.field_names = inp.validate_and_convert('field', fn, 'list(string)')

        # Get the functions and verify the function spaces
        for fn in self.field_names:
            func_3d = simulation.get_data(fn)
            V = func_3d.function_space()
            fam = V.ufl_element().family()
            deg = V.ufl_element().degree()
            if self.family is None:
                self.family = fam
                self.degree = deg
            elif fam != self.family or deg != self.degree:
                ocellaris_error(
                    'Mismatching function spaces in PlainProbe %s' % self.name,
                    'All functions must have the same function space. '
                    + '%s is %r but %r was expected' % (fn, (fam, deg), (self.family, self.degree)),
                )

        # Create the slice
        self.slice = FunctionSlice(self.plane_point, self.plane_normal, V, xlim, ylim, zlim)
        prefix = simulation.input.get_value('output/prefix', '', 'string')

        # Get the XDMF file name (also ensures it does not exist)
        fn = '%s_slice_%s.xdmf' % (prefix, self.name)
        self.file_name = get_xdmf_file_name(simulation, fn)

        if simulation.rank == 0:
            V_2d = self.slice.slice_function_space
            mesh_2d = V_2d.mesh()
            simulation.log.info('        Created 2D mesh with %r cells' % mesh_2d.num_cells())
            simulation.log.info('        Creating XDMF file %s' % self.file_name)
            self.xdmf_file = dolfin.XDMFFile(dolfin.MPI.comm_self, self.file_name)
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

        # Add field to list of IO plotters
        inp_key = probe_input.basepath + 'write_interval'
        simulation.io.add_plotter(self.write_field, inp_key, WRITE_INTERVAL)

    def write_field(self):
        """
        Find and output the plane probe
        """
        for fn, func_2d in zip(self.field_names, self.funcs_2d):
            func_3d = self.simulation.get_data(fn)
            self.slice.get_slice(func_3d, func_2d)
            if self.simulation.rank == 0:
                self.xdmf_file.write(func_2d, self.simulation.time)


class FunctionSlice:
    def __init__(self, pt, n, V3d, xlim=None, ylim=None, zlim=None):
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
        mesh_2d, cell_origins = make_cut_plane_mesh(pt, n, V3d.mesh(), xlim, ylim, zlim)

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


def make_cut_plane_mesh(pt, n, mesh3d, xlim=None, ylim=None, zlim=None):
    """
    Returns a 2D mesh of the intersection of a 3D mesh and a plane. The result
    can optionally be restricted to the intersection of the plane and a cube
    with sides parallel to the axes by specifying one or more of the coordinate
    limits as tuples, they are by default None which means _lim = (-inf, inf).

    * pt: a point in the plane
    * n: a normal vector to the plane. Does not need to be a unit normal
    * mesh3d: the 3D mesh to be intersected by the plane
    * xlim, ylim, zlim: each of these must be a tuple of two numbers or None

    This function assumes that the 3D mesh consists solely of tetrahedra which
    gives a 2D mesh of triangles
    """
    # Get results on this rank and send to root process
    rank_results = get_points_in_plane(pt, n, mesh3d)
    rank_results = split_cells(rank_results)
    rank_results = limit_plane(rank_results, xlim, ylim, zlim)
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
        for cell_id, subcells in res.items():
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


def limit_plane(unlimited, xlim=None, ylim=None, zlim=None, eps=1e-8):
    """
    Given a set of triangles in a 2D plane, limit the extents of the plane to a
    given quadratic area which sides are parallel to the axes
    """
    if xlim is None and ylim is None and zlim is None:
        return unlimited

    # Remove / alter cells that are outside the coordinate limits
    results = OrderedDict()
    for cell_id, subcells in unlimited.items():
        subcells = limit_triangles(subcells, xlim, ylim, zlim, eps)
        if subcells:
            results[cell_id] = subcells

    return results


def limit_triangles(triangles, xlim, ylim, zlim, eps=1e-8):
    """
    Given a list of triangles, return a new list of triangles which covers the
    same area, limited to the given axis limits (i.e. the area of the triangles
    returned is always equal to or lower than the area of the input triangles.
    An ampty list can be returned if no area exists inside the axes limits
    """
    for d, lim in enumerate((xlim, ylim, zlim)):
        if lim is None:
            continue

        # Check lower limit
        new_triangles = []
        for cell_coords in triangles:
            below = [1 if c[d] < lim[0] - eps else 0 for c in cell_coords]
            num_below = sum(below)

            if num_below == 0:
                # Cell completely inside the bounds
                new_triangles.append(cell_coords)
                continue

            elif num_below in (1, 2):
                # Cell partially inside the bounds
                if num_below == 1:
                    # Get the index of the one vertex that is below lim
                    i = 0 if below[0] == 1 else (1 if below[1] == 1 else 2)
                else:
                    # Get the index of the one vertex that is not below lim
                    i = 0 if below[0] == 0 else (1 if below[1] == 0 else 2)

                # Coordinates of the vertices
                c0 = cell_coords[(i + 0) % 3]
                c1 = cell_coords[(i + 1) % 3]
                c2 = cell_coords[(i + 2) % 3]

                # Coordinates of the crossing points
                f01 = (lim[0] - c0[d]) / (c1[d] - c0[d])
                f02 = (lim[0] - c0[d]) / (c2[d] - c0[d])
                c01 = (
                    c0[0] * (1 - f01) + c1[0] * f01,
                    c0[1] * (1 - f01) + c1[1] * f01,
                    c0[2] * (1 - f01) + c1[2] * f01,
                )
                c02 = (
                    c0[0] * (1 - f02) + c2[0] * f02,
                    c0[1] * (1 - f02) + c2[1] * f02,
                    c0[2] * (1 - f02) + c2[2] * f02,
                )

                # Create new triangles that are inside the bounds
                if num_below == 1:
                    # Split into two new triangles
                    new_triangles.append((c1, c01, c2))
                    new_triangles.append((c2, c01, c02))
                else:
                    new_triangles.append((c0, c01, c02))

            elif num_below == 3:
                # Cell completely outside the bounds
                continue
        triangles = new_triangles

        # Check upper limit
        new_triangles = []
        for cell_coords in triangles:
            above = [1 if c[d] > lim[1] + eps else 0 for c in cell_coords]
            num_above = sum(above)

            if num_above == 0:
                # Cell completely inside the bounds
                new_triangles.append(cell_coords)
                continue

            elif num_above in (1, 2):
                # Cell partially inside the bounds
                if num_above == 1:
                    # Get the index of the one vertex that is above lim
                    i = 0 if above[0] == 1 else (1 if above[1] == 1 else 2)
                else:
                    # Get the index of the one vertex that is not above lim
                    i = 0 if above[0] == 0 else (1 if above[1] == 0 else 2)

                # Coordinates of the vertices
                c0 = cell_coords[(i + 0) % 3]
                c1 = cell_coords[(i + 1) % 3]
                c2 = cell_coords[(i + 2) % 3]

                # Coordinates of the crossing points
                f01 = (lim[1] - c0[d]) / (c1[d] - c0[d])
                f02 = (lim[1] - c0[d]) / (c2[d] - c0[d])
                c01 = (
                    c0[0] * (1 - f01) + c1[0] * f01,
                    c0[1] * (1 - f01) + c1[1] * f01,
                    c0[2] * (1 - f01) + c1[2] * f01,
                )
                c02 = (
                    c0[0] * (1 - f02) + c2[0] * f02,
                    c0[1] * (1 - f02) + c2[1] * f02,
                    c0[2] * (1 - f02) + c2[2] * f02,
                )

                # Create new triangles that are inside the bounds
                if num_above == 1:
                    # Split into two new triangles
                    new_triangles.append((c1, c01, c2))
                    new_triangles.append((c2, c01, c02))
                else:
                    new_triangles.append((c0, c01, c02))

            elif num_above == 3:
                # Cell completely outside the bounds
                continue
        triangles = new_triangles

    def has_area(cell_coords):
        c0, c1, c2 = cell_coords
        u = (c1[0] - c0[0], c1[1] - c0[1], c1[2] - c0[2])
        v = (c2[0] - c0[0], c2[1] - c0[1], c2[2] - c0[2])
        # Cross product to find area of trapezoid, squared
        areaish = (
            (u[1] * v[2] - u[2] * v[1]) ** 2
            + (u[2] * v[0] - u[0] * v[2]) ** 2
            + (u[0] * v[1] - u[1] * v[0]) ** 2
        )
        return areaish > eps ** 2

    return [cell_coords for cell_coords in triangles if has_area(cell_coords)]


def split_cells(unsplit):
    """
    Split non-triangles into triangles
    """
    results = OrderedDict()
    for cell_id, cell_coords in unsplit.items():
        N = len(cell_coords)
        if N == 4:
            results[cell_id] = [cell_coords[1:], cell_coords[:-1]]
        elif N == 3:
            results[cell_id] = [cell_coords]
        else:
            raise NotImplementedError('Expected elements with 3 or 4 ' 'vertices, not %d' % N)
    return results


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

    def sign(x):
        return -1 if x < 0 else 1

    n = plane_coefficients[:3]  # normal to the plane

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

    return w / (w.dot(w)) ** 0.5


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
