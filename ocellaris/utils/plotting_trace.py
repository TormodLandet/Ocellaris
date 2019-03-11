# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
import numpy


def plot_matplotlib_dgt(func, **kwargs):
    """
    Plot scalar discontinuous Galerkin Trace functions in 2D
    """
    # Get information about the underlying function space
    function_space = func.function_space()
    family = func.ufl_element().family()
    mesh = function_space.mesh()
    ndim = mesh.geometry().dim()

    # Check that the function is of a supported type
    if not family == 'HDiv Trace':
        raise NotImplementedError('Trace function plotter does not support %r' % family)
    if not ndim == 2:
        raise NotImplementedError(
            'Trace function plotter only supports 2D, you gave %dD' % ndim
        )

    # Get dof, area and midpoint for each facet
    facet_dofmap, areas, midpoints = facet_info_and_dofmap(function_space)
    avg_facet_length = numpy.average(areas)
    Nfacet = len(facet_dofmap)

    # Build matplotlib plolygon data
    coords = mesh.coordinates()
    verts = []
    for cell in dolfin.cells(mesh):
        corners = cell.entities(0)
        xs = [coords[i, 0] for i in corners]
        ys = [coords[i, 1] for i in corners]
        verts.append(list(zip(xs, ys)))

    # Rearange function values according to the dofmap
    scalars = numpy.zeros(Nfacet, float)
    funcvals = func.vector().get_local()

    for i, facet in enumerate(dolfin.facets(mesh)):
        fidx = facet.index()
        scalars[i] = funcvals[facet_dofmap[fidx]]

    # Size of scatterpoints and plot axes (we need to set axis limits for this plot)
    radius = kwargs.pop('radius', avg_facet_length / 8)
    xlim = kwargs.pop('xlim', (coords[:, 0].min(), coords[:, 0].max()))
    ylim = kwargs.pop('ylim', (coords[:, 1].min(), coords[:, 1].max()))

    scalar_locs = numpy.array(midpoints)
    return plot_triangulation_scatter(
        verts, scalar_locs, scalars, radius, xlim, ylim, **kwargs
    )


def plot_triangulation_scatter(vertices, points, scalars, radius, xlim, ylim, **kwargs):
    """
    Plot the 2D geometry on top of a scatterplot with scalars
    """
    from matplotlib.collections import PolyCollection
    from matplotlib import pyplot as plt

    ax = plt.gca()

    # Find  the length of the scalar circle radius in pixels
    # We need to have the correct xlim and ylim first!
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ptx, pty = ax.transData.transform(
        [(radius, 0), (0, radius)]
    ) - ax.transData.transform((0, 0))
    rad_px = min(ptx[0], pty[1])
    area_px = numpy.pi * rad_px ** 2

    vmin = scalars.min()
    vmax = scalars.max()
    if vmin == vmax == 0:
        vmax = 1

    # Plot the scalars with colors
    patches = ax.scatter(
        points[:, 0],
        points[:, 1],
        s=area_px,
        c=scalars,
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        **kwargs
    )

    # Add grid data above the scalars
    polys = PolyCollection(
        vertices, facecolors=(1, 1, 1, 0), edgecolors='black', linewidths=0.25
    )
    ax.add_collection(polys)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    return patches


def facet_info_and_dofmap(V):
    """
    Return three lists

      1) A list which maps facet index to dof
      2) A list which maps facet index to facet area
      2) A list which maps facet index to facet midpoint
    """
    mesh = V.mesh()
    dofmap = V.dofmap()

    ndim = V.ufl_cell().topological_dimension()

    # Loop through cells and get dofs for each cell
    facet_dofs = [None] * mesh.num_facets()
    facet_area = [None] * mesh.num_facets()
    facet_midp = [None] * mesh.num_facets()
    for cell in dolfin.cells(mesh):
        dofs = dofmap.cell_dofs(cell.index())
        facet_idxs = cell.entities(ndim - 1)

        # Only works for functions with one dof on each facet
        assert len(dofs) == len(facet_idxs)

        # Loop through connected facets and store dofs for each facet
        for i, fidx in enumerate(facet_idxs):
            facet_dofs[fidx] = dofs[i]
            facet_area[fidx] = cell.facet_area(i)
            mp = dolfin.Facet(mesh, fidx).midpoint()
            facet_midp[fidx] = (mp.x(), mp.y())

    return facet_dofs, facet_area, facet_midp
