# encoding: utf-8
from __future__ import division
import sys
from xml.etree.ElementTree import Element, SubElement, register_namespace, tostring


def read_iso_surface_file(file_name):
    with open(file_name, 'rt') as f:
        # Read header
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()

        field_name = line1[31:-7]
        value = float(line2[9:])
        dim = int(line3[8:])
        
        assert dim == 2
        times = []
        lines = []
        
        tline = f.readline()
        while tline:
            wt = tline.split()
            time = float(wt[1])
            nlines = int(wt[3])
            
            tlines = []
            for _ in range(nlines):
                xvals = [float(v) for v in f.readline().split()]
                yvals = [float(v) for v in f.readline().split()]
                zvals = [float(v) for v in f.readline().split()]
                tlines.append((xvals, yvals, zvals))
            
            times.append(time)
            lines.append(tlines)
            tline = f.readline()
        
        return field_name, value, dim, times, lines
    

class X3DMesh(object):
    def __init__(self):
        self.connectivity = []
        self.coordinates = []
    
    def _write(self, scene_element):
        group = SubElement(scene_element, 'Group')
        shape = SubElement(group, 'Shape')
        ifs = SubElement(shape, 'IndexedFaceSet')
        coords = SubElement(ifs, 'Coordinate')
        
        ifs.set('coordIndex', ' -1 '.join(' '.join(str(v) for v in vs) for vs in self.connectivity) + ' -1')
        coords.set('point', ' '.join(' '.join(str(x) for x in xs) for xs in self.coordinates))


class X3DFile(object):
    def __init__(self):
        self.meshes = []
    
    def write(self, x3d_file_name):
        # Doctype and root element with namespaces etc
        XSD = 'http://www.w3.org/2001/XMLSchema-instance'
        register_namespace('xsd', XSD)
        x3d = Element('X3D', {'profile': 'Interchange', 'version': '3.3'})
        x3d.set('{%s}noNamespaceSchemaLocation' % XSD, 'http://www.web3d.org/specifications/x3d-3.3.xsd')
        
        # Metadata
        head = SubElement(x3d, 'head')
        SubElement(head, 'meta', {'name': 'generator', 'content': 'Ocellaris iso_line_to_xd3.py'})
        
        # The Scene
        scene = SubElement(x3d, 'Scene')
        #SubElement(scene, 'Background', {'skyColor': '1 1 1'})
        #SubElement(scene, 'Viewpoint', {'orientation': '0 -1 0 0.53', 'position': '-2.28 0.29 4.06'})
        
        for mesh in self.meshes:
            mesh._write(scene)
    
        xml = '<?xml version="1.0" encoding="UTF-8"?>\n'
        xml += '<!DOCTYPE X3D PUBLIC "ISO//Web3D//DTD X3D 3.3//EN" "http://www.web3d.org/specifications/x3d-3.3.dtd">\n'
        xml += tostring(x3d, encoding='utf-8')
        
        with open(x3d_file_name, 'wt') as out:
            out.write(xml)


def iso_line_to_x3d(iso_file_name, x3d_file_name):
    field_name, value, dim, times, lines = read_iso_surface_file(iso_file_name)
    print 'Processing field %s with iso-line at %g, dim = %r' % (field_name, value, dim)
    
    x3d = X3DFile()
    
    for itime, t in enumerate(times):
        tlines = lines[itime]
        x3d.meshes = []
        
        for xs, ys, zs in tlines:
            mesh = X3DMesh()
            x3d.meshes.append(mesh)
            
            # Add vertex coordinates
            mesh.coordinates += [(x, 0, y) for x, y, _z in zip(xs, ys, zs)]
            mesh.coordinates += [(x, 1, y) for x, y, _z in zip(xs, ys, zs)]
            
            # Construct connectivity
            N = len(xs)
            for i in range(N-1):
                mesh.connectivity.append((i, i + 1, N + i + 1, N + i))
        
        fn = x3d_file_name % itime
        print 'Time %10.5f: %s' % (t, fn)
        x3d.write(fn)


if __name__ == '__main__':
    iso_line_file_name = sys.argv[1]
    x3d_file_name = sys.argv[2]
    iso_line_to_x3d(iso_line_file_name, x3d_file_name)
