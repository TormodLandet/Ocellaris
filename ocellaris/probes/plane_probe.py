import dolfin
import numpy
from ocellaris.utils import init_mesh_geometry


def make_cut_plane_mesh(pt, n, mesh3d):
    """
    Returns a 2D mesh of the intersection of a 3D mesh and a plane
    
    * pt: a point in the plane
    * n: a normal vector to the plane. Does not need to be a unit normal
    * mesh3d: the 3D mesh to be intersected by the plane
    
    This function assumes that the 3D mesh consists solely of tetrahedra which
    gives a 2D mesh of triangles 
    """
    res = get_points_in_plane(pt, n, mesh3d)
    comm = mesh3d.mpi_comm()
    
    results = comm.rank, res
    all_results = comm.gather(results)
    
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
    
    results = {}
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
                cell_points.append((pos[0], pos[1], pos[2], cid))
        
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
