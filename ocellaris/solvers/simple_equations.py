# encoding: utf8
from __future__ import division
from collections import deque
import numpy
from petsc4py import PETSc
import dolfin
from dolfin import dx, div, grad, dot, jump, avg, dS, Constant
from ocellaris.utils import timeit
from ocellaris.solver_parts import define_penalty


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
        
        # Storage for forms
        self.eqAs, self.eqBs, self.eqCs, self.eqDs, self.eqE = [], [], [], [], 0
        
        # Create UFL forms
        self.define_coupled_equation()
        
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
        
    def calculate_penalties(self, nu, no_coeff=False):
        """
        Calculate SIPG penalty
        """
        mpm = self.simulation.multi_phase_model
        mesh = self.simulation.data['mesh']
        
        if no_coeff:
            mu_min = mu_max = 1.0
        else:
            mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
        
        P = self.simulation.data['Vu'].ufl_element().degree()
        penalty_dS = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0)
        penalty_ds = penalty_dS*2
        self.simulation.log.info('    DG SIP penalty:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
        
        if False and self.velocity_continuity_factor_D12 is not None:
            D12 = Constant([self.velocity_continuity_factor_D12]*self.simulation.ndim)
        else:
            D12 = Constant([0, 0])
        
        if False and self.pressure_continuity_factor != 0:
            h = self.simulation.data['h']
            h = Constant(1.0)
            D11 = avg(h/nu)*Constant(self.pressure_continuity_factor)
        else:
            D11 = None
        
        return Constant(penalty_dS), Constant(penalty_ds), D11, D12
    
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mpm = sim.multi_phase_model
        mesh = sim.data['mesh']
        u_conv = sim.data['u_conv']
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        
        # The trial and test functions
        ndim = self.simulation.ndim
        ulist = [dolfin.TrialFunction(Vu) for _ in range(ndim)]
        vlist = [dolfin.TestFunction(Vu) for _ in range(ndim)]
        u = dolfin.as_vector(ulist)
        v = dolfin.as_vector(vlist)
        p = dolfin.TrialFunction(Vp)
        q = dolfin.TestFunction(Vp)
        
        c1, c2, c3 = sim.data['time_coeffs']
        dt = sim.data['dt']
        g = sim.data['g']
        n = dolfin.FacetNormal(mesh)
        
        # Fluid properties
        rho = mpm.get_density(0)
        nu = mpm.get_laminar_kinematic_viscosity(0)
        mu = mpm.get_laminar_dynamic_viscosity(0)
        
        # Hydrostatic pressure correction
        if self.include_hydrostatic_pressure:
            p += sim.data['p_hydrostatic']
        
        # Penalties
        penalty_dS, penalty_ds, D11, D12 = self.calculate_penalties(nu)
            
        # Upwind and downwind velocities
        w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2.0
        w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2.0
        
        # Momentum equations
        for d in range(sim.ndim):
            up = sim.data['up%d' % d]
            upp = sim.data['upp%d' % d]
            eqA = eqB = eqC = 0
            
            # Divergence free criterion
            # ∇⋅u = 0
            if self.incompressibility_flux_type == 'central':
                u_hat_p = avg(u[d])
            elif self.incompressibility_flux_type == 'upwind':
                assert self.use_grad_q_form, 'Upwind only implemented for grad_q_form'
                switch = dolfin.conditional(dolfin.gt(w_nU('+'), 0.0), 1.0, 0.0)
                u_hat_p = switch*u[d]('+') + (1 - switch)*u[d]('-')
            
            if self.use_grad_q_form:
                eqC -= u[d]*q.dx(d)*dx
                eqC += (u_hat_p + D12[d]*jump(u, n))*jump(q)*n[d]('+')*dS
            else:
                eqC += q*u[d].dx(d)*dx
                eqC -= (avg(q) - dot(D12, jump(q, n)))*jump(u[d])*n[d]('+')*dS
            
            # Time derivative
            # ∂(ρu)/∂t
            eqA += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
            
            # Convection:
            # -w⋅∇(ρu)
            flux_nU = u[d]*w_nU
            flux = jump(flux_nU)
            eqA -= u[d]*dot(grad(rho*v[d]), u_conv)*dx
            eqA += flux*jump(rho*v[d])*dS
            
            # Stabilizing term when w is not divergence free
            eqA += 1/2*div(u_conv)*u[d]*v[d]*dx
            
            # Diffusion:
            # -∇⋅μ∇u
            eqA += mu*dot(grad(u[d]), grad(v[d]))*dx
            
            # Symmetric Interior Penalty method for -∇⋅μ∇u
            eqA -= avg(mu)*dot(n('+'), avg(grad(u[d])))*jump(v[d])*dS
            eqA -= avg(mu)*dot(n('+'), avg(grad(v[d])))*jump(u[d])*dS
            
            # Symmetric Interior Penalty coercivity term
            eqA += penalty_dS*jump(u[d])*jump(v[d])*dS
            
            # Pressure
            # ∇p
            if self.use_grad_p_form:
                eqB += v[d]*p.dx(d)*dx
                eqB -= (avg(v[d]) + D12[d]*jump(v, n))*jump(p)*n[d]('+')*dS
            else:
                eqB -= p*v[d].dx(d)*dx
                eqB += (avg(p) - dot(D12, jump(p, n)))*jump(v[d])*n[d]('+')*dS
            
            # Pressure continuity stabilization. Needed for equal order discretization
            if D11 is not None:
                eqB += D11*dot(jump(p, n), jump(q, n))*dS
            
            # Body force (gravity)
            # ρ g
            eqA -= rho*g[d]*v[d]*dx
            
            # Other sources
            for f in sim.data['momentum_sources']:
                eqA -= f[d]*v[d]*dx
            
            # Dirichlet boundary
            dirichlet_bcs = sim.data['dirichlet_bcs'].get('u%d' % d, [])
            for dbc in dirichlet_bcs:
                u_bc = dbc.func()
                
                # Divergence free criterion
                if self.use_grad_q_form:
                    eqC += q*u_bc*n[d]*dbc.ds()
                else:
                    eqC -= q*u[d]*n[d]*dbc.ds()
                    eqC += q*u_bc*n[d]*dbc.ds()
                
                # Convection
                eqA += rho*u[d]*w_nU*v[d]*dbc.ds()
                eqA += rho*u_bc*w_nD*v[d]*dbc.ds()
                
                # SIPG for -∇⋅μ∇u
                eqA -= mu*dot(n, grad(u[d]))*v[d]*dbc.ds()
                eqA -= mu*dot(n, grad(v[d]))*u[d]*dbc.ds()
                eqA += mu*dot(n, grad(v[d]))*u_bc*dbc.ds()
                
                # Weak Dirichlet
                eqA += penalty_ds*(u[d] - u_bc)*v[d]*dbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    eqB += p*v[d]*n[d]*dbc.ds()
            
            # Neumann boundary
            neumann_bcs = sim.data['neumann_bcs'].get('u%d' % d, [])
            for nbc in neumann_bcs:
                # Divergence free criterion
                if self.use_grad_q_form:
                    eqC += q*u[d]*n[d]*nbc.ds()
                else:
                    eqC -= q*u[d]*n[d]*nbc.ds()
                
                # Convection
                eqA += rho*u[d]*w_nU*v[d]*nbc.ds()
                
                # Diffusion
                eqA -= mu*nbc.func()*v[d]*nbc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    eqB += p*v[d]*n[d]*nbc.ds()
            
            # Outlet boundary
            for obc in sim.data['outlet_bcs']:
                # Divergence free criterion
                if self.use_grad_q_form:
                    eqC += q*u[d]*n[d]*obc.ds()
                else:
                    eqC -= q*u[d]*n[d]*obc.ds()
                
                # Convection
                eqA += rho*u[d]*w_nU*v[d]*obc.ds()
                
                # Diffusion
                mu_dudn = p*n[d]
                eqA -= mu_dudn*v[d]*obc.ds()
                
                # Pressure
                if not self.use_grad_p_form:
                    p_ = mu*dot(dot(grad(u), n), n)
                    eqB += p_*n[d]*v[d]*obc.ds()
            
            eqA, eqD1 = dolfin.system(eqA)
            eqB, eqD2 = dolfin.system(eqB)
            eqC, eqEd = dolfin.system(eqC) 
            self.eqAs.append(eqA)
            self.eqBs.append(eqB)
            self.eqCs.append(eqC)
            self.eqDs.append(eqD1 + eqD2)
            self.eqE += eqEd
    
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
            At = create_block_matrix(self.simulation, Vu, 1)
            Ati = create_block_matrix(self.simulation, Vu, 1)
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
            At = create_block_matrix(self.simulation, Vu, self.block_partitions)
            Ati = create_block_matrix(self.simulation, Vu, self.block_partitions)
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


def create_block_matrix(simulation, V, blocks=1):
    """
    Create a sparse PETSc matrix to hold dense blocks that are larger than
    the normal DG block diagonal mass matrices (super-cell dense blocks)
    Based on code from Miro: https://fenicsproject.org/qa/11647
    """
    mesh = simulation.data['mesh']
    dm = V.dofmap()
    
    if isinstance(blocks, int):
        nnz_max = blocks
    else:
        # Create block diagonal matrix
        nnz_max = 0
        for _cells, dofs, _dof_idx in blocks:
            nnz_max = max(nnz_max, len(dofs))
    
    comm = mesh.mpi_comm().tompi4py() 
    mat = PETSc.Mat()
    mat.create(comm)
    
    # Set local and global size of the matrix
    sizes = [dm.index_map().size(dolfin.IndexMap.MapSize_OWNED), 
             dm.index_map().size(dolfin.IndexMap.MapSize_GLOBAL)]
    mat.setSizes([sizes, sizes])
    
    # Set matrix type
    mat.setType('aij')
    mat.setPreallocationNNZ(nnz_max)
    mat.setUp()
    
    # Map from local rows to global rows
    lgmap = map(int, dm.tabulate_local_to_global_dofs())
    lgmap = PETSc.LGMap().create(lgmap, comm=comm)
    mat.setLGMap(lgmap, lgmap)
    
    return dolfin.PETScMatrix(mat)


EQUATION_SUBTYPES = {
    'Default': SimpleEquations,
}