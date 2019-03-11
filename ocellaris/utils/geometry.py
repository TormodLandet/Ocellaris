# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
from collections import namedtuple


# Compact way to store the information, may turn into classes later if needed
CellInfo = namedtuple('CellData', 'volume midpoint')
FacetInfo = namedtuple('FacetData', 'area midpoint normal on_boundary')


def init_connectivity(simulation):
    """
    Initialize the needed connectivity data
    """
    mesh = simulation.data['mesh']
    ndim = simulation.ndim

    if ndim == 2:
        # In 2D: cell = face, facet = edge

        # Connectivity from vertex to face and vice versa
        mesh.init(0, 2)
        mesh.init(2, 0)

        # Connectivity from face to edge and vice versa
        mesh.init(2, 1)
        mesh.init(1, 2)

        # Connectivity from vertex to edge and vice versa
        mesh.init(0, 1)
        mesh.init(1, 0)

        simulation.data['connectivity_VC'] = mesh.topology()(0, 2)
        simulation.data['connectivity_CV'] = mesh.topology()(2, 0)
        simulation.data['connectivity_FC'] = mesh.topology()(1, 2)
        simulation.data['connectivity_CF'] = mesh.topology()(2, 1)
        simulation.data['connectivity_FV'] = mesh.topology()(1, 0)
        simulation.data['connectivity_VF'] = mesh.topology()(0, 1)

    else:
        # Connectivity from vertex to cell and vice versa
        mesh.init(0, 3)
        mesh.init(3, 0)

        # Connectivity from cell to face and vice versa
        mesh.init(3, 2)
        mesh.init(2, 3)

        simulation.data['connectivity_VC'] = mesh.topology()(0, 3)
        simulation.data['connectivity_CV'] = mesh.topology()(3, 0)
        simulation.data['connectivity_FC'] = mesh.topology()(2, 3)
        simulation.data['connectivity_CF'] = mesh.topology()(3, 2)
        simulation.data['connectivity_FV'] = mesh.topology()(2, 0)
        simulation.data['connectivity_VF'] = mesh.topology()(0, 2)


def precompute_cell_data(simulation):
    """
    Get cell volume and midpoint in an easy to use format
    """
    mesh = simulation.data['mesh']
    ndim = simulation.ndim

    cell_info = [None] * mesh.num_cells()
    for cell in dolfin.cells(mesh, 'all'):
        mp = cell.midpoint()
        if ndim == 2:
            midpoint = numpy.array([mp.x(), mp.y()], float)
        else:
            midpoint = numpy.array([mp.x(), mp.y(), mp.z()], float)

        volume = cell.volume()
        cell_info[cell.index()] = CellInfo(volume, midpoint)

    simulation.data['cell_info'] = cell_info


def precompute_facet_data(simulation):
    """
    Get facet normal and areas in an easy to use format
    """
    mesh = simulation.data['mesh']
    conFC = simulation.data['connectivity_FC']
    ndim = simulation.ndim
    cell_info = simulation.data['cell_info']

    # Get the facet areas from the cells
    areas = {}
    for cell in dolfin.cells(mesh, 'all'):
        # Get the connected facets
        facet_idxs = cell.entities(ndim - 1)

        # Loop over connected facets and get the area
        for i, fidx in enumerate(facet_idxs):
            a = cell.facet_area(i)
            if fidx in areas:
                assert a == areas[fidx]
            else:
                areas[fidx] = a

    # Loop over facets and gather the required information
    facet_info = [None] * mesh.num_facets()
    for facet in dolfin.facets(mesh, 'all'):
        fidx = facet.index()

        mp = facet.midpoint()
        if ndim == 2:
            midpoint = numpy.array([mp.x(), mp.y()], float)
        else:
            midpoint = numpy.array([mp.x(), mp.y(), mp.z()], float)

        # Find one cell connected to this facet. There can be one or two
        # connected cells, we only need the first one
        connected_cells = conFC(fidx)
        icell0 = connected_cells[0]
        on_boundary = len(connected_cells) == 1

        # Midpoint of local cell 0
        cell0_mp = cell_info[icell0].midpoint

        # Vector from local cell midpoint to face midpoint
        vec0 = midpoint - cell0_mp

        # Find a normal pointing out from cell 0
        normalpt = facet.normal()
        if ndim == 2:
            normal = numpy.array([normalpt.x(), normalpt.y()], float)
        else:
            normal = numpy.array([normalpt.x(), normalpt.y(), normalpt.z()], float)
        if numpy.dot(vec0, normal) < 0:
            normal *= -1

        area = areas[fidx]
        facet_info[fidx] = FacetInfo(area, midpoint, normal, on_boundary)

    simulation.data['facet_info'] = facet_info
