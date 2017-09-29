from collections import deque
import numpy
import dolfin
from ocellaris.utils import timeit, split_form_into_matrix, create_block_matrix
from .coupled_equations import define_dg_equations


class SimpleEquations(object):
    def __init__(self, simulation, use_stress_divergence_form, use_grad_p_form,
                 use_grad_q_form, use_lagrange_multiplicator, 
                 include_hydrostatic_pressure, incompressibility_flux_type):
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
        self.num_elements_in_block = 0
        self.lump_diagonal = False
        self.block_partitions = None
        
        assert self.incompressibility_flux_type in ('central', 'upwind')
        
        # Discontinuous or continuous elements
        Vu_family = simulation.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # We do not currently support all possible options
        assert self.vel_is_discontinuous
        assert not self.simulation.mesh_morpher.active
        assert not self.use_lagrange_multiplicator
        assert not self.use_stress_divergence_form
        
        # Create UFL forms
        self.eqAs, self.eqBs, self.eqCs, self.eqDs, self.eqE = [], [], [], [], 0
        self.define_simple_equations()
        
        # Storage for assembled matrices
        self.As = [None] * simulation.ndim
        self.A_tildes = [None] * simulation.ndim
        self.A_tilde_invs = [None] * simulation.ndim
        self.Bs = [None] * simulation.ndim
        self.Cs = [None] * simulation.ndim
        self.Ds = [None] * simulation.ndim
        self.E = None
        self.E_star = None
        self.L = [None] * simulation.ndim
    
    def define_simple_equations(self):
        """
        Setup weak forms for SIMPLE form
        """
        sim = self.simulation
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        
        # The trial and test functions in a coupled space (to be split)
        func_spaces = [Vu] * sim.ndim + [Vp]
        e_mixed = dolfin.MixedElement([fs.ufl_element() for fs in func_spaces])
        Vcoupled = dolfin.FunctionSpace(sim.data['mesh'], e_mixed)
        tests = dolfin.TestFunctions(Vcoupled)
        trials = dolfin.TrialFunctions(Vcoupled)
        
        # Split into components
        v = dolfin.as_vector(tests[:-1])
        u = dolfin.as_vector(trials[:-1])
        q = tests[-1]
        p = trials[-1]
        lm_trial = lm_test = None
        
        # Define the full coupled form and split it into subforms depending
        # on the test and trial functions
        eq = define_dg_equations(u, v, p, q, lm_trial, lm_test, self.simulation,
                                 include_hydrostatic_pressure=self.include_hydrostatic_pressure,
                                 incompressibility_flux_type=self.incompressibility_flux_type,
                                 use_grad_q_form=self.use_grad_q_form,
                                 use_grad_p_form=self.use_grad_p_form,
                                 use_stress_divergence_form=self.use_stress_divergence_form)
        mat, vec = split_form_into_matrix(eq, Vcoupled, Vcoupled)
        
        # There is no p*q form, this is a saddle point system
        assert mat[-1,-1] is None, 'Found p-q coupling, this is not a saddle point system!'
        
        # Store the forms
        for i in range(sim.ndim):
            self.eqAs.append(mat[i,i]) 
            self.eqBs.append(mat[i,-1])
            self.eqCs.append(mat[-1,i])
            self.eqDs.append(vec[i])
            self.eqE = vec[-1]
            
            # Check that off diagonal terms are zero
            for j in range(sim.ndim):
                if i == j:
                    continue
                assert mat[i,j] is None, 'No coupling between velocity components supported in SIMPLE solver!'
    
    @timeit
    def assemble_matrices(self, d, reassemble=False):
        # Equations, matrices and flag to indicate reassembly needed
        eqs_and_matrices = ((self.eqAs, self.As, True),
                            (self.eqBs, self.Bs, reassemble),
                            (self.eqCs, self.Cs, reassemble))
        
        # Assemble A, B and C matrices
        for eqs, Ms, reas in eqs_and_matrices:
            if Ms[d] is None:
                Ms[d] = dolfin.as_backend_type(dolfin.assemble(eqs[d]))
            elif reas:
                dolfin.assemble(eqs[d], tensor=Ms[d])
        
        # Assemble Ã and Ã_inv matrices
        Nelem = self.num_elements_in_block
        if Nelem == 0:
            At, Ati = self.assemble_A_tilde_diagonal(d)
        elif Nelem == 1:
            At, Ati = self.assemble_A_tilde_single_element(d)
        else:
            At, Ati = self.assemble_A_tilde_multi_element(d, Nelem)
        At.apply('insert'); Ati.apply('insert')
        self.A_tildes[d], self.A_tilde_invs[d] = At, Ati
        
        return (self.As[d], self.A_tildes[d], self.A_tilde_invs[d],
                self.Bs[d], self.Cs[d])
    
    @timeit
    def assemble_A_tilde_diagonal(self, d):
        """
        Assemble diagonal Ã and Ã_inv matrices
        """
        Vu = self.simulation.data['Vu']
        Aglobal = dolfin.as_backend_type(self.As[d])
        
        if self.A_tildes[d] is None:
            At = create_block_matrix(Vu, 1)
            Ati = create_block_matrix(Vu, 1)
            self.u_diag = dolfin.Vector(self.simulation.data['u0'].vector())
        else:
            At = self.A_tildes[d]
            Ati = self.A_tilde_invs[d]
        At.zero()
        Ati.zero()
        
        if self.lump_diagonal:
            self.u_diag[:] = 1.0
            self.u_diag.apply('insert')
            self.u_diag = Aglobal * self.u_diag
        else:
            Aglobal.get_diagonal(self.u_diag)
        
        At.set_diagonal(self.u_diag)
        inv_diag = 1.0/self.u_diag.get_local()
        self.u_diag.set_local(inv_diag)
        Ati.set_diagonal(self.u_diag)
        
        return At, Ati
    
    @timeit
    def assemble_A_tilde_single_element(self, d):
        """
        Assemble block diagonal Ã and Ã_inv matrices where the blocks
        are the dofs in a single element
        """
        Aglobal = dolfin.as_backend_type(self.As[d])
        if self.A_tildes[d] is None:
            At = dolfin.PETScMatrix(Aglobal)
            Ati = dolfin.PETScMatrix(Aglobal)
        else:
            At = self.A_tildes[d]
            Ati = self.A_tilde_invs[d]
        At.zero()
        Ati.zero()
        
        dm = self.simulation.data['Vu'].dofmap()
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
    def assemble_A_tilde_multi_element(self, d, Nelem):
        """
        Assemble block diagonal Ã and Ã_inv matrices where the blocks
        are the dofs of N elememts a single element
        """
        Vu = self.simulation.data['Vu']
        if self.block_partitions is None:
            self.block_partitions = create_block_partitions(self.simulation, Vu, Nelem)
            self.simulation.log.info('SIMPLE solver with %d cell blocks found %d blocks in total'
                                     % (Nelem, len(self.block_partitions)))
        
        Aglobal = dolfin.as_backend_type(self.As[d])
        if self.A_tildes[d] is None:
            #At = dolfin.PETScMatrix(Aglobal)
            At = create_block_matrix(Vu, self.block_partitions)
            Ati = create_block_matrix(Vu, self.block_partitions)
        else:
            At = self.A_tildes[d]
            Ati = self.A_tilde_invs[d]
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
    def assemble_D(self, d):
        if self.Ds[d] is None:
            self.Ds[d] = dolfin.assemble(self.eqDs[d])
        else:
            dolfin.assemble(self.eqDs[d], tensor=self.Ds[d])
        return self.Ds[d]

    @timeit
    def assemble_E(self):
        if self.E is None:
            self.E = dolfin.assemble(self.eqE)
        else:
            dolfin.assemble(self.eqE, tensor=self.E)
        return self.E
    
    @timeit
    def assemble_E_star(self, u_star):
        if self.E_star is None:
            self.E_star = dolfin.Vector(self.simulation.data['p'].vector())
        E_star = self.E_star
        E_star.zero()
        
        # Divergence of u*, C⋅u*
        for d in range(self.simulation.ndim):
            E_star.axpy(1.0, self.Cs[d]*u_star[d].vector())
        
        # Subtract the original RHS of C⋅u = e
        E_star.axpy(-1.0, self.assemble_E())
        
        return E_star


def create_block_partitions(simulation, V, Ncells):
    """
    Create super-cell partitions of Ncells cells each 
    """
    mesh = simulation.data['mesh']
    
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

    # Get dofs per cell
    dm = V.dofmap()
    Ndof = dm.cell_dofs(0).shape[0]
    N = Ncells*Ndof
    
    # Partition all local cells into super-cells 
    picked = [False]*num_cells_owned
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


EQUATION_SUBTYPES = {
    'Default': SimpleEquations,
}