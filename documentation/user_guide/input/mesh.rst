.. _inp_mesh:

Mesh
====

You can specify simple geometries using FEniCS DOLFIN built in mesh generators,
and also load a mesh from file. For realistic cases using something like gmsh
to generate meshes is recommended. The meshio_ program can be used to convert
between different mesh file formats and also loading these formats directly,
see below.

.. _meshio: https://github.com/nschloe/meshio


Simple geometries
-----------------

Example: 2D rectangle

.. code-block:: yaml
        
    mesh:
        type: Rectangle
        Nx: 64
        Ny: 64
        diagonal: left/right  # defaults to 'right'
        startx: 0             # defaults to 0
        endx:   2             # defaults to 1
        # you can also give starty and endy

Example: 3D box

.. code-block:: yaml
        
    mesh:
        type: Box
        Nx: 64
        Ny: 64
        Nz: 15
        startx: 0  # defaults to 0
        endx:   2  # defaults to 1
        # you can also give starty and endy, startz and endz

Example: 2D disc

.. code-block:: yaml
        
    mesh:
        type: UnitDisc
        N: 20
        degree: 1  # defaults to 1 (degree of mesh elements)


Mesh file formats
-----------------

**Example:** using meshio_ to load all its supported formats (RECOMMENDED)

.. code-block:: yaml
        
    mesh:
        type: meshio
        mesh_file: mesh.msh
        meshio_type: gmsh

The supported formats (as of November 2018) can be found `in this list 
<https://github.com/nschloe/meshio/blob/8289814be4f714b6d6000e173ab6697d1f35655f/meshio/helpers.py#L130>`_
in the meshio source on github.

**Example:** legacy DOLFIN XML format

.. code-block:: yaml
        
    mesh:
        type: XML
        mesh_file: mesh.xml
        facet_region_file: regions.xml  # not required

Ocellaris will look for the xml files first as absolute paths, then as paths
relative to the current working directory and last as paths relative to the
directory of the input file. If it cannot find the file in any of these
places you will get an error message and Ocellaris will quit.

A sample mesh xml file and facet marker file is included in the ``demo/files``
directory. The mesh ``ocellaris_mesh.xml.gz`` and the facet regions
``ocellaris_facet_regions.xml.gz``. You can load these files without unzipping
them. The *flow around Ocellaris* demo shows how it is done.

**Example:** XDMF format

.. code-block:: yaml
        
    mesh:
        type: XDMF
        mesh_file: mesh.xdmf

**Example:** Ocellaris HDF5 restart file format

.. code-block:: yaml
        
    mesh:
        type: HDF5
        mesh_file: ocellaris_savepoint000010.h5

This will only load the mesh and (possibly) facet regions. You can also start
the simulation from a restart file instead of an input file. Then the mesh *and*
the function values from that save point are used, allowing you to restart the
simulation more or less like it was never stopped.


Moving the mesh
---------------

Ocellaris can move the mesh right after it has been created or read from file.
To move the mesh in order to refine, skew, scale, rotate or translate it you
must specify a C++ description of the mesh *displacement* from the initial
position (which was specified in the input file or in the loaded mesh file).

An example is the following 140 meter long 2D wave tank which is 10 m high. To
refine the mesh in the y-direction such that it is finest around ``x[1] = 7``
meters—where the free surface is to be located—a function is specified which
is zero on the boundaries (to avoid changing the domain size) and non-zero in
the interior in order to move the nodes closer to the free surface. No refinement
is performed in the x-direction (``x[0]``).

.. code-block:: yaml
        
    mesh:
        type: Rectangle
        Nx: 140
        Ny: 20
        endx: 140
        endy: 20
        move: ['0', '0.0297619048*pow(x[1], 3) - 0.520833333*pow(x[1], 2) + 2.23214286*x[1] + 3.55271368e-15']

In order to develop and check the mesh refinement function it can be beneficial
to generate and plot it, e.g., using matplotlib in jupyter or using similar
interactive tools. The above refinement was developed using polynomial fitting
in numpy::

    from matplotlib import pyplot
    import numpy
    
    # Find a polynomial that refines the mesh
    y_target = [0, 4, 7.5, 10]
    dy_target = [0, 2.5, 0, 0]  # zero at the boundary
    P = numpy.polyfit(y_target, dy_target, 3)
    
    # Realise the polynomial
    y = numpy.linspace(0, 10, 20)
    dy = numpy.polyval(P, y)
    
    # Plot the results
    for ypos in (y + dy):
        pyplot.plot([0, 1], [ypos, ypos], '-k', lw=1)'
    pyplot.axhline(7, c='b', ls=':')
    pyplot.axhline(6, c='b', ls=':', lw=1)
    pyplot.axhline(8, c='b', ls=':', lw=1)

For more complicated meshes it is recommended to perform mesh grading and other
mesh operation in an external mesh generator such as gmsh. 
There is also some (not much used, hence possibly buggy) support for ALE where
the mesh moves every timestep, but that is not covered by the ``mesh`` section
of the input file.