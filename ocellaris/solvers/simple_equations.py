# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from collections import deque
import numpy
import dolfin
from ocellaris.utils import timeit, split_form_into_matrix, create_block_matrix
from .coupled_equations import define_dg_equations


class SimpleEquations(object):
    def __init__(
        self,
        simulation,
        use_stress_divergence_form,
        use_grad_p_form,
        use_grad_q_form,
        use_lagrange_multiplicator,
        include_hydrostatic_pressure,
        incompressibility_flux_type,
        num_elements_in_block,
        lump_diagonal,
        a_tilde_is_mass=False,
    ):
        """
        This class assembles the coupled Navier-Stokes equations as a set of
        matrices and vectors

            | A  B |   | u |   | D |
            |      | . |   | = |   |
            | C  0 |   | p |   | 0 |

        There is also the vector E = C u, since this will not be zero until
        the iterations have converged. In addition we have Ã and Ãinv which
        are approximations to the A and Ainv matrices in such a way that the
        inverse is easy to compute. We use Ã = I * (1 + time derivative part)

        :type simulation: ocellaris.Simulation
        """
        self.simulation = simulation
        self.use_stress_divergence_form = use_stress_divergence_form
        self.use_grad_p_form = use_grad_p_form
        self.use_grad_q_form = use_grad_q_form
        self.use_lagrange_multiplicator = use_lagrange_multiplicator
        self.include_hydrostatic_pressure = include_hydrostatic_pressure
        self.incompressibility_flux_type = incompressibility_flux_type
        self.num_elements_in_block = num_elements_in_block
        self.lump_diagonal = lump_diagonal
        self.block_partitions = None
        self.a_tilde_is_mass = a_tilde_is_mass

        # Show some configuration info
        log = simulation.log.info
        nelem = self.num_elements_in_block
        if nelem == 0:
            log('    Using diagonal Ã matrix')
        elif nelem == 1:
            log('    Using block-diagonal Ã matrix')
        else:
            log('    Using %d-element block-diagonal Ã matrix' % nelem)
        if self.a_tilde_is_mass:
            log('    Using only mass matrix when constructing Ã')
        else:
            log('    Using full A matrix when constructing Ã')
        if self.lump_diagonal:
            log('    Lumping diagonal constructing Ã')

        # We do not currently support all possible options
        assert self.incompressibility_flux_type in ('central', 'upwind')
        assert not self.simulation.mesh_morpher.active
        assert not self.use_lagrange_multiplicator
        assert not self.use_stress_divergence_form

        # Create UFL forms
        self.eqA = self.eqB = self.eqC = self.eqD = self.eqE = None
        self.define_simple_equations()

        # Storage for assembled matrices
        self.A = self.A_tilde = self.A_tilde_inv = None
        self.B = self.C = self.D = self.E = None
        self.L = None

    def define_simple_equations(self):
        """
        Setup weak forms for SIMPLE form
        """
        sim = self.simulation
        self.Vuvw = sim.data['uvw_star'].function_space()
        Vp = sim.data['Vp']

        # The trial and test functions in a coupled space (to be split)
        func_spaces = [self.Vuvw, Vp]
        e_mixed = dolfin.MixedElement([fs.ufl_element() for fs in func_spaces])
        Vcoupled = dolfin.FunctionSpace(sim.data['mesh'], e_mixed)
        tests = dolfin.TestFunctions(Vcoupled)
        trials = dolfin.TrialFunctions(Vcoupled)

        # Split into components
        v = dolfin.as_vector(tests[0][:])
        u = dolfin.as_vector(trials[0][:])
        q = tests[-1]
        p = trials[-1]
        lm_trial = lm_test = None

        # Define the full coupled form and split it into subforms depending
        # on the test and trial functions
        eq = define_dg_equations(
            u,
            v,
            p,
            q,
            lm_trial,
            lm_test,
            self.simulation,
            include_hydrostatic_pressure=self.include_hydrostatic_pressure,
            incompressibility_flux_type=self.incompressibility_flux_type,
            use_grad_q_form=self.use_grad_q_form,
            use_grad_p_form=self.use_grad_p_form,
            use_stress_divergence_form=self.use_stress_divergence_form,
        )
        mat, vec = split_form_into_matrix(eq, Vcoupled, Vcoupled, check_zeros=True)

        # Check matrix and vector shapes and that the matrix is a saddle point matrix
        assert mat.shape == (2, 2)
        assert vec.shape == (2,)
        assert mat[-1, -1] is None, 'Found p-q coupling, this is not a saddle point system!'

        # Store the forms
        self.eqA = mat[0, 0]
        self.eqB = mat[0, 1]
        self.eqC = mat[1, 0]
        self.eqD = vec[0]
        self.eqE = vec[1]

        if self.eqE is None:
            self.eqE = dolfin.TrialFunction(Vp) * dolfin.Constant(0) * dolfin.dx

        if self.a_tilde_is_mass:
            # The mass matrix. Consistent with the implementation in define_dg_equations
            rho = sim.multi_phase_model.get_density(0)
            c1 = sim.data['time_coeffs'][0]
            dt = sim.data['dt']
            eqM = rho * c1 / dt * dolfin.dot(u, v) * dolfin.dx
            matM, _vecM = split_form_into_matrix(eqM, Vcoupled, Vcoupled, check_zeros=True)
            self.eqM = dolfin.Form(matM[0, 0])
            self.M = None

    @timeit
    def assemble_matrices(self, reassemble=False):
        # Assemble self.A, self.B and self.C matrices from the
        # weak forms self.eqA, self.eqB and self.eqC
        for name, reas in [('A', True), ('B', reassemble), ('C', reassemble)]:
            M = getattr(self, name)
            eq = getattr(self, 'eq' + name)
            if M is None:
                M = dolfin.as_backend_type(dolfin.assemble(eq))
                setattr(self, name, M)
            elif reas:
                dolfin.assemble(eq, tensor=M)

        if self.a_tilde_is_mass:
            if self.M is None:
                self.M = dolfin.as_backend_type(dolfin.assemble(self.eqM))
            else:
                dolfin.assemble(self.eqM, tensor=self.M)

        # Assemble Ã and Ã_inv matrices
        Nelem = self.num_elements_in_block
        if Nelem == 0:
            At, Ati = self.assemble_A_tilde_diagonal()
        elif Nelem == 1:
            At, Ati = self.assemble_A_tilde_single_element()
        else:
            At, Ati = self.assemble_A_tilde_multi_element(Nelem)
        At.apply('insert')
        Ati.apply('insert')
        self.A_tilde, self.A_tilde_inv = At, Ati

        return self.A, self.A_tilde, self.A_tilde_inv, self.B, self.C

    @timeit
    def assemble_A_tilde_diagonal(self):
        """
        Assemble diagonal Ã and Ã_inv matrices
        """
        uvw = self.simulation.data['uvw_star']
        Vuvw = uvw.function_space()
        Aglobal = self.M if self.a_tilde_is_mass else self.A
        if self.A_tilde is None:
            At = create_block_matrix(Vuvw, 'diag')
            Ati = create_block_matrix(Vuvw, 'diag')
            self.u_diag = dolfin.Vector(uvw.vector())
        else:
            At = self.A_tilde
            Ati = self.A_tilde_inv
        At.zero()
        Ati.zero()

        if self.lump_diagonal:
            self.u_diag[:] = 1.0
            self.u_diag.apply('insert')
            self.u_diag = Aglobal * self.u_diag
        else:
            Aglobal.get_diagonal(self.u_diag)

        # Setup the diagonal A matrix
        At.set_diagonal(self.u_diag)

        # Setup the inverse of the diagonal A matrix
        inv_diag = 1.0 / self.u_diag.get_local()
        self.u_diag.set_local(inv_diag)
        self.u_diag.apply('insert')
        Ati.set_diagonal(self.u_diag)

        return At, Ati

    @timeit
    def assemble_A_tilde_single_element(self):
        """
        Assemble block diagonal Ã and Ã_inv matrices where the blocks
        are the dofs in a single element
        """
        Aglobal = self.M if self.a_tilde_is_mass else self.A
        if self.A_tilde is None:
            At = Aglobal.copy()
            Ati = Aglobal.copy()
        else:
            At = self.A_tilde
            Ati = self.A_tilde_inv
        At.zero()
        Ati.zero()

        dm = self.Vuvw.dofmap()
        N = dm.cell_dofs(0).shape[0]
        Alocal = numpy.zeros((N, N), float)

        # Loop over cells and get the block diagonal parts (should be moved to C++)
        for cell in dolfin.cells(self.simulation.data['mesh'], 'regular'):
            # Get global dofs
            istart = Aglobal.local_range(0)[0]
            dofs = dm.cell_dofs(cell.index()) + istart

            # Get block diagonal part of A, invert and insert into approximations
            Aglobal.get(Alocal, dofs, dofs)
            Alocal_inv = numpy.linalg.inv(Alocal)
            At.set(Alocal, dofs, dofs)
            Ati.set(Alocal_inv, dofs, dofs)
        return At, Ati

    @timeit
    def assemble_A_tilde_multi_element(self, Nelem):
        """
        Assemble block diagonal Ã and Ã_inv matrices where the blocks
        are the dofs of N elememts a single element
        """
        if self.block_partitions is None:
            self.block_partitions = create_block_partitions(self.simulation, self.Vuvw, Nelem)
            self.simulation.log.info(
                'SIMPLE solver with %d cell blocks found %d blocks in total'
                % (Nelem, len(self.block_partitions))
            )

        Aglobal = self.M if self.a_tilde_is_mass else self.A
        if self.A_tilde is None:
            block_dofs = [dofs for _, dofs, _ in self.block_partitions]
            At = create_block_matrix(self.Vuvw, block_dofs)
            Ati = At.copy()
        else:
            At = self.A_tilde
            Ati = self.A_tilde_inv
        At.zero()
        Ati.zero()

        # Loop over super-cells and get the block diagonal parts (should be moved to C++)
        istart = Aglobal.local_range(0)[0]
        for _cells, dofs, _dof_idx in self.block_partitions:
            global_dofs = dofs + istart
            N = len(dofs)
            Ablock = numpy.zeros((N, N), float)
            Aglobal.get(Ablock, global_dofs, global_dofs)
            Ablock_inv = numpy.linalg.inv(Ablock)

            At.set(Ablock, dofs, global_dofs)
            Ati.set(Ablock_inv, dofs, global_dofs)
        return At, Ati

    @timeit
    def assemble_D(self):
        if self.D is None:
            self.D = dolfin.assemble(self.eqD)
        else:
            dolfin.assemble(self.eqD, tensor=self.D)
        return self.D

    @timeit
    def assemble_E(self):
        if self.E is None:
            self.E = dolfin.assemble(self.eqE)
        else:
            dolfin.assemble(self.eqE, tensor=self.E)
        return self.E

    @timeit
    def assemble_E_star(self, u_star):
        # Divergence of u*, C⋅u*
        E_star = self.C * u_star.vector()

        # Subtract the original RHS of C⋅u = e
        E_star.axpy(-1.0, self.assemble_E())

        return E_star


def create_block_partitions(simulation, V, Ncells):
    """
    Create super-cell partitions of Ncells cells each
    """
    mesh = simulation.data['mesh']
    dm = V.dofmap()

    # Construct a cell connectivity mapping
    con_CF = simulation.data['connectivity_CF']
    con_FC = simulation.data['connectivity_FC']
    con_CFC = {}
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    for icell in range(num_cells_owned):
        for ifacet in con_CF(icell):
            for inb in con_FC(ifacet):
                if inb != icell:
                    con_CFC.setdefault(icell, []).append(inb)

    # Partition all local cells into super-cells
    picked = [False] * num_cells_owned
    partitions = []
    for icell in range(num_cells_owned):
        # Make sure the cell is not part of an existing supercell
        if picked[icell]:
            continue

        # Find candidate cells to join this supercell and
        # extend the candidate set by breadth first search
        super_cell = [icell]
        picked[icell] = True
        candidates = deque(con_CFC[icell])
        while candidates and len(super_cell) < Ncells:
            icand = candidates.popleft()
            if picked[icand]:
                continue

            super_cell.append(icand)
            picked[icand] = True
            candidates.extend(con_CFC[icand])

        # Get the dofs of our super-cell
        # Will contain duplicates if V is not DG
        dofs = []
        for isel in super_cell:
            dofs.extend(dm.cell_dofs(isel))

        # Map dofs to indices in local block matrix
        dof_idx = {}
        for i, dof in enumerate(dofs):
            dof_idx[dof] = i

        dofs = numpy.array(dofs, numpy.intc)
        partitions.append((super_cell, dofs, dof_idx))

    return partitions


EQUATION_SUBTYPES = {'Default': SimpleEquations, 'DG': SimpleEquations}
