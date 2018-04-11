import dolfin
import numpy
from ocellaris.probes.plane_probe import (make_cut_plane_mesh, get_points_in_plane,
                                          get_plane_normal, get_point_sides,
                                          get_plane_coefficients)


def test_cut_mesh():
    mesh = dolfin.UnitCubeMesh(2, 2, 2)
    mesh.init(3, 0)
    pt = (0.1, 0.1, 0.1)
    n = (0, 0, 1)
    mesh2d, cell_origins = make_cut_plane_mesh(pt, n, mesh)
    
    import pprint; pprint.pprint(cell_origins)
    
    if False and mesh2d is not None:
        import matplotlib; matplotlib.use('Agg')
        from matplotlib import pyplot
        fig = pyplot.figure()
        dolfin.plot(mesh2d)
        fig.savefig('test_cut_mesh.png')


def test_points_in_mesh():
    mesh = dolfin.UnitCubeMesh(1, 1, 1)
    mesh.init(3, 0)
    
    for i, pt, n in [(0, (0.5, 0.5, 0.5), (3, 0, 0)),
                     (1, (0.5, 0.5, 0.5), (0, 1, 0)),
                     (2, (0.5, 0.5, 0.5), (0, 0, 9))]:
        results = get_points_in_plane(pt, n, mesh)
        import pprint; pprint.pprint(results)
        
        for _cid, pts in results.items():
            for pt in pts:
                assert pt[i] == 0.5
                assert pt[(i + 1) % 3] in (0, 0.5, 1)
                assert pt[(i + 2) % 3] in (0, 0.5, 1)


def test_plane_coefficients():
    for pts in [(( 0,  0,  0), ( 1,  0,  0), ( 0,  0,  1)),
                ((-2,  1,  1), ( 0,  2,  3), ( 1,  0, -1))]:
        p1, p2, p3 = pts
        n = get_plane_normal(p1, p2, p3)
        c = get_plane_coefficients(p1, n)
        
        f = numpy.min([abs(v) for v in c if abs(v) > 1e-8])
        allpts = []
        print('\npoints', pts, '\nnormal', n, '\nans', c * 1 / f)
        
        # Verify plane coefficients by checking that the points are coinciding
        # with the computed plane
        for p in (p1, p2, p3):
            s = numpy.dot(p, c[:3]) + c[3]
            assert abs(s) < 1e-15
            allpts.append(p)
        
        # Other points should not be incident to the computed plane
        for p in [(10,  -10,  10), n + p1, n - p2, n + p3]:
            assert abs(numpy.dot(p, c[:3]) + c[3]) > 1e-6
            allpts.append(p)
        
        allpts = numpy.array(allpts, float)
        sides = get_point_sides(c, allpts)
        s2 = abs(sides) < 1e-8
        assert (s2 == [True, True, True, False, False, False, False]).all()
