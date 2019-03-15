# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
from dolfin import Timer, Constant, Cell, cells


class SlopeLimiterBoundaryConditions(object):
    BC_TYPE_NOT_ON_BOUNDARY = 0
    BC_TYPE_DIRICHLET = 1
    BC_TYPE_NEUMANN = 2
    BC_TYPE_ROBIN = 3
    BC_TYPE_OTHER = 4

    def __init__(self, simulation, field_name, dof_region_marks, V):
        """
        This class helps slope limiting of cells adjacent to the boundary by
        providing either values of the function and its derivatives at the
        boundary

        @param dict dof_region_marks: As created by get_dof_region_marks()
        """
        assert field_name is not None
        self.simulation = simulation

        # Velocity BCs can be given for u itself (coupled components like slip
        # BCs) or for u0, u1 and u2 separately (e.g., for pure Dirichlet BCs)
        if field_name[-1] in '012':
            self.field_names = [field_name[:-1], field_name]
        else:
            self.field_names = [field_name]

        self.function_space = V
        self.active = False
        self.set_dof_region_marks(dof_region_marks)
        self._dof_to_cell = None
        self._allready_warned = set()

    def set_dof_region_marks(self, dof_region_marks):
        """
        In case boundary regions move around the marks can be updated here

        @param dict dof_region_marks: map from dof to list of regions
            containing the dof (multiple regions happens in corners etc)
        """

        self.dof_region_marks = {}  # Map from dof to ONE region
        self.region_dofs = {}  # Map from region number to list of dofs

        for dof, regions in dof_region_marks.items():
            region = regions[-1]  # in case of multiple regions pick the last
            self.dof_region_marks[dof] = region
            self.region_dofs.setdefault(region, []).append(dof)

    def activate(self, active=True):
        self.active = active

    def _warn(self, message):
        if message not in self._allready_warned:
            self.simulation.log.warning(message)
            self._allready_warned.add(message)

    def get_bcs(self, weak_dof_values=None):
        """
        Get bc type and bc value for each dof in the function space
        If a weak_dof_values array is given then this is used for Dirichlet BCs
        instead of the strong BCs given by the BC function (these will typically
        differ by a small bit)
        """
        sim = self.simulation
        im = self.function_space.dofmap().index_map()
        num_owned_dofs = im.size(im.MapSize.OWNED)
        boundary_dof_type = numpy.zeros(num_owned_dofs, numpy.intc)
        boundary_dof_value = numpy.zeros(num_owned_dofs, float)

        if not self.active:
            return boundary_dof_type, boundary_dof_value

        # This is potentially slow, so we time this code
        timer = Timer("Ocellaris get slope limiter boundary conditions")

        # Collect BCs - field name for u0 can be u (FreeSlip) and u0 (CppCodedValue)
        dirichlet = {}
        neumann = {}
        robin = {}
        slip = {}
        for field_name in self.field_names:
            # Collect Dirichlet BCs for this field
            for bc in sim.data['dirichlet_bcs'].get(field_name, []):
                region_number = bc.subdomain_id - 1
                dirichlet[region_number] = bc

            # Collect Neumann BCs for this field
            for bc in sim.data['neumann_bcs'].get(field_name, []):
                region_number = bc.subdomain_id - 1
                neumann[region_number] = bc

            # Collect Robin BCs for this field
            for bc in sim.data['robin_bcs'].get(field_name, []):
                region_number = bc.subdomain_id - 1
                robin[region_number] = bc

            # Collect Slip BCs for this field
            for bc in sim.data['slip_bcs'].get(field_name, []):
                region_number = bc.subdomain_id - 1
                slip[region_number] = bc

        fname = ', '.join(self.field_names)
        regions = sim.data['boundary']
        for region_number, dofs in self.region_dofs.items():
            boundary_region = regions[region_number]

            # Get the BC object
            if region_number in dirichlet:
                bc_type = self.BC_TYPE_DIRICHLET
                value = dirichlet[region_number].func()
            elif region_number in neumann:
                bc_type = self.BC_TYPE_NEUMANN
                value = neumann[region_number].func()
            elif region_number in robin:
                bc_type = self.BC_TYPE_ROBIN
                value = robin[region_number].dfunc()
            elif region_number in slip:
                value = None
                for dof in dofs:
                    boundary_dof_type[dof] = self.BC_TYPE_OTHER
                continue
            else:
                self._warn(
                    'WARNING: Slope limiter found no BC for field %s '
                    'in region %s' % (fname, boundary_region.name)
                )
                continue

            if bc_type == self.BC_TYPE_DIRICHLET and weak_dof_values is not None:
                # Take values from a field which has (weak) Dirichlet BCs applied
                for dof in dofs:
                    boundary_dof_type[dof] = bc_type
                    boundary_dof_value[dof] = weak_dof_values[dof]

            elif isinstance(value, Constant):
                # Get constant value
                val = value.values()
                assert val.size == 1
                val = val[0]

                for dof in dofs:
                    boundary_dof_type[dof] = bc_type
                    boundary_dof_value[dof] = val

            elif hasattr(value, 'eval_cell'):
                # Get values from an Expression of some sort
                dof_to_cell = self._get_dof_to_cell_mapping()
                mesh = self.function_space.mesh()
                val = numpy.zeros(1, float)
                for dof in dofs:
                    cid, coords = dof_to_cell[dof]
                    cell = Cell(mesh, cid)
                    value.eval_cell(val, coords, cell)
                    boundary_dof_type[dof] = bc_type
                    boundary_dof_value[dof] = val[0]

            else:
                self._warn(
                    'WARNING: Field %s has unsupported limiter BC %r in '
                    'region %s' % (fname, type(value), boundary_region.name)
                )

        timer.stop()
        return boundary_dof_type, boundary_dof_value

    def _get_dof_to_cell_mapping(self):
        """
        Mapping from dof to containing cell and coordinates of location
        """
        if self._dof_to_cell is not None:
            return self._dof_to_cell

        V = self.function_space
        mesh = V.mesh()
        gdim = mesh.geometry().dim()
        dofs_x = V.tabulate_dof_coordinates().reshape((-1, gdim))
        dm = V.dofmap()

        self._dof_to_cell = dof_to_cell = [None] * dofs_x.shape[0]
        for cell in cells(mesh):
            cid = cell.index()
            for dof in dm.cell_dofs(cid):
                dof_coords = numpy.array(dofs_x[dof], float)
                dof_to_cell[dof] = (cid, dof_coords)

        return dof_to_cell
