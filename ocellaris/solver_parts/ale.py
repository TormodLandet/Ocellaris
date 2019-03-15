# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import timeit, ocellaris_interpolate


class MeshMorpher(object):
    def __init__(self, simulation):
        """
        Class to handle prescribed or runtime evaluated mesh morphing
        """
        self.simulation = simulation
        self.active = False

        # The user can give a mesh velocity function to simulate a piston or similar
        prescribed_velocity_input = simulation.input.get_value('mesh/prescribed_velocity', None)
        if prescribed_velocity_input is not None:
            self.setup_prescribed_velocity(prescribed_velocity_input)

        # Fields that will be interpolated to the new location
        self.previous_fields = [
            'up0',
            'up1',
            'up2',
            'upp0',
            'upp1',
            'upp2',
            'uppp0',
            'uppp1',
            'uppp2',
        ]

    def setup(self):
        """
        Create mesh velocity and deformation functions
        """
        sim = self.simulation
        assert self.active is False, 'Trying to setup mesh morphing twice in the same simulation'

        # Store previous cell volumes
        mesh = sim.data['mesh']
        Vcvol = dolfin.FunctionSpace(mesh, 'DG', 0)
        sim.data['cvolp'] = dolfin.Function(Vcvol)

        # The function spaces for mesh velocities and displacements
        Vmesh = dolfin.FunctionSpace(mesh, 'CG', 1)
        Vmesh_vec = dolfin.VectorFunctionSpace(mesh, 'CG', 1)
        sim.data['Vmesh'] = Vmesh

        # Create mesh velocity functions
        u_mesh = []
        for d in range(sim.ndim):
            umi = dolfin.Function(Vmesh)
            sim.data['u_mesh%d' % d] = umi
            u_mesh.append(umi)
        u_mesh = dolfin.as_vector(u_mesh)
        sim.data['u_mesh'] = u_mesh

        # Create mesh displacement vector function
        self.displacement = dolfin.Function(Vmesh_vec)
        self.assigners = [dolfin.FunctionAssigner(Vmesh_vec.sub(d), Vmesh) for d in range(sim.ndim)]
        self.active = True

    def setup_prescribed_velocity(self, input_dict):
        """
        The mesh nodes will be updated onece every time step
        """
        # Read and verify input
        assert input_dict['type'] == 'CppCodedValue'
        self.cpp_codes = input_dict['cpp_code']
        assert isinstance(self.cpp_codes, list) and len(self.cpp_codes) == self.simulation.ndim

        # Setup mesh morphing every time step
        self.setup()
        self.simulation.hooks.add_pre_timestep_hook(
            lambda timestep_number, t, dt: self.update_prescribed_mesh_velocity(),
            'MeshMorpher - update prescribed mesh velocity',
        )

    @timeit
    def update_prescribed_mesh_velocity(self):
        """
        Move the mesh according to prescribed velocities
        """
        sim = self.simulation
        Vmesh = sim.data['Vmesh']

        for d, cpp_code in enumerate(self.cpp_codes):
            description = 'initial conditions for mesh vel %d' % d
            func = sim.data['u_mesh%d' % d]

            # Update the mesh velocity functions
            ocellaris_interpolate(sim, cpp_code, description, Vmesh, func)

        self.morph_mesh()

    def morph_mesh(self):
        """
        Move the mesh and update cached geometry information. It is assumed
        that the mesh velocities u_mesh0, u_mesh1 etc are already populated
        """
        sim = self.simulation

        # Get the mesh displacement
        for d in range(sim.ndim):
            umi = sim.data['u_mesh%d' % d]
            self.assigners[d].assign(self.displacement.sub(d), umi)
        self.displacement.vector()[:] *= sim.dt

        # Save the cell volumes before morphing
        mesh = sim.data['mesh']
        cvolp = sim.data['cvolp']
        dofmap_cvol = cvolp.function_space().dofmap()
        for cell in dolfin.cells(mesh):
            dofs = dofmap_cvol.cell_dofs(cell.index())
            assert len(dofs) == 1
            cvolp.vector()[dofs[0]] = cell.volume()

        # Move the mesh according to the given displacements
        dolfin.ALE.move(mesh, self.displacement)
        mesh.bounding_box_tree().build(mesh)
        sim.update_mesh_data(connectivity_changed=False)
