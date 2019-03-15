# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os
import dolfin


# Default values, can be changed in the input file
XDMF_FLUSH = True


class XDMFFileIO:
    def __init__(self, simulation):
        """
        Xdmf output using dolfin.XDMFFile
        """
        self.simulation = simulation
        self.extra_functions = []
        self.xdmf_file = None

    def close(self):
        """
        Close any open xdmf file
        """
        if self.xdmf_file is not None:
            self.xdmf_file.close()
        self.xdmf_file = None

    def write(self):
        """
        Write a file that can be used for visualization. The fluid fields will
        be automatically downgraded (interpolated) into something that
        dolfin.XDMFFile can write, typically linear CG elements.
        """
        with dolfin.Timer('Ocellaris save xdmf'):
            if self.xdmf_file is None:
                self._setup_xdmf()
            self._write_xdmf()
        return self.file_name

    def _setup_xdmf(self):
        """
        Create XDMF file object
        """
        sim = self.simulation
        xdmf_flush = sim.input.get_value('output/xdmf_flush', XDMF_FLUSH, 'bool')

        fn = sim.input.get_output_file_path('output/xdmf_file_name', '.xdmf')
        file_name = get_xdmf_file_name(sim, fn)

        sim.log.info('    Creating XDMF file %s' % file_name)
        comm = sim.data['mesh'].mpi_comm()
        self.file_name = file_name
        self.xdmf_file = dolfin.XDMFFile(comm, file_name)
        self.xdmf_file.parameters['flush_output'] = xdmf_flush
        self.xdmf_file.parameters['rewrite_function_mesh'] = False
        self.xdmf_file.parameters['functions_share_mesh'] = True
        self.xdmf_first_output = True

        def create_vec_func(V):
            "Create a vector function from the components"
            family = V.ufl_element().family()
            degree = V.ufl_element().degree()
            cd = sim.data['constrained_domain']
            V_vec = dolfin.VectorFunctionSpace(
                sim.data['mesh'], family, degree, constrained_domain=cd
            )
            vec_func = dolfin.Function(V_vec)
            assigner = dolfin.FunctionAssigner(V_vec, [V] * sim.ndim)
            return vec_func, assigner

        # XDMF cannot save functions given as "as_vector(list)"
        self._vel_func, self._vel_func_assigner = create_vec_func(sim.data['Vu'])
        self._vel_func.rename('u', 'Velocity')
        if sim.mesh_morpher.active:
            self.xdmf_file.parameters['rewrite_function_mesh'] = True
            self._mesh_vel_func, self._mesh_vel_func_assigner = create_vec_func(sim.data['Vmesh'])
            self._mesh_vel_func.rename('u_mesh', 'Velocity of the mesh')

    def _write_xdmf(self):
        """
        Write plot files for Paraview and similar applications
        """
        sim = self.simulation
        t = float(sim.time)

        if self.xdmf_first_output:
            bm = sim.data['boundary_marker']
            self.xdmf_file.write(bm)

        # Write the fluid velocities
        vel = sim.data.get('up', sim.data['u'])
        self._vel_func_assigner.assign(self._vel_func, list(vel))
        self.xdmf_file.write(self._vel_func, t)

        # Write the mesh velocities (used in ALE calculations)
        if sim.mesh_morpher.active:
            self._mesh_vel_func_assigner.assign(self._mesh_vel_func, list(sim.data['u_mesh']))
            self.xdmf_file.write(self._mesh_vel_func, t)

        # Write scalar functions
        for name in ('p', 'p_hydrostatic', 'c', 'rho'):
            if name in sim.data:
                func = sim.data[name]
                if isinstance(func, dolfin.Function):
                    self.xdmf_file.write(func, t)

        # Write extra functions
        for func in self.extra_functions:
            self.xdmf_file.write(func, t)

        self.xdmf_first_output = False


def get_xdmf_file_name(simulation, file_name_suggestion):
    """
    Deletes any previous files with the same name unless
    the simulation is restarted. If will then return a new
    name to avoid overwriting the existing files
    """
    base = os.path.splitext(file_name_suggestion)[0]
    if simulation.restarted:
        base = base + '_restarted_%08d' % simulation.timestep

    file_name = base + '.xdmf'
    file_name2 = base + '.h5'

    # Remove any existing files with the same base file name
    if os.path.isfile(file_name) and simulation.rank == 0:
        simulation.log.info('    Removing existing XDMF file %s' % file_name)
        os.remove(file_name)
    if os.path.isfile(file_name2) and simulation.rank == 0:
        simulation.log.info('    Removing existing XDMF file %s' % file_name2)
        os.remove(file_name2)

    return file_name
