from __future__ import division
import dolfin
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error
from . import Solver, register_solver, define_advection_problem, define_poisson_problem

@register_solver('IPCS')
class SolverIPCS(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the pressure-velocity splitting
        scheme IPCS (Incremental Pressure Correction Scheme)
        """
        self.simulation = sim = simulation
        mesh = sim.data['mesh']
        
        # Test for PETSc
        if not dolfin.has_linear_algebra_backend("PETSc"):
            report_error('Missing PETSc',
                         'DOLFIN has not been configured with PETSc '
                         'which is needed by Ocellaris.')
        dolfin.parameters["linear_algebra_backend"] = "PETSc"
        
        # Function spaces
        Vu = simulation.data['Vu']
        Vp = simulation.data['Vp']
        
        # Create velocity functions. Keep both component and vector forms
        uvec, upvec, uppvec, u_conv, u_star = [], [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            sim.data['u_star%d' % d] = us = dolfin.Function(Vu)
            uvec.append(u)
            upvec.append(up)
            uppvec.append(upp)
            u_conv.append(uc)
            u_star.append(us)
        sim.data['u'] = dolfin.as_vector(uvec)
        sim.data['up'] = dolfin.as_vector(upvec)
        sim.data['upp'] = dolfin.as_vector(uppvec)
        sim.data['u_conv'] = u_conv = dolfin.as_vector(u_conv)
        sim.data['u_star'] = u_star = dolfin.as_vector(u_star)
        
        # Create pressure function
        sim.data['p'] = p = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)
        
        # Mesh parameters
        n = dolfin.FacetNormal(mesh)
        #h = dolfin.CellSize(mesh)
        
        # Physical properties
        rho = sim.data['rho']
        nu = sim.data['nu']
        g = sim.data['g']
        self.dt = dolfin.Constant(1.0)
        
        # Get convection schemes for the velocity
        conv_schemes = []
        for d in range(sim.ndim):
            conv_scheme_name = sim.input['convection']['u']['convection_scheme']
            conv_scheme = get_convection_scheme(conv_scheme_name)(simulation, 'u%d' % d)
            conv_schemes.append(conv_scheme)
        self.convection_schemes = conv_schemes
        
        # Define the momentum prediction equations
        h = 1/mesh.num_cells()
        penalty = 1.0/h
        self.eqs_mom_pred = []
        for d in range(sim.ndim):
            beta = conv_schemes[d].blending_function
            f = -1/rho*p.dx(d) + g[d]
            a1, L1 = define_advection_problem(Vu, upvec[d], uppvec[d], 
                                              u_conv, n, beta, self.dt)
            a2, L2 = define_poisson_problem(Vu, -nu, f, n, penalty, [])
            eq = a1+a2, L1+L2
            self.eqs_mom_pred.append(eq)
        
        # Define the pressure correction equation
        penalty = 1.0/h
        f = 3/(2*self.dt)*dolfin.nabla_div(u_star)
        neumann_bcs = self.simulation.data['neumann_bcs'].get('p', [])
        self.eq_pressure = define_poisson_problem(Vp, 1/rho, f, n, penalty, neumann_bcs)
        
    def update_convection(self, t, dt):
        """
        Update terms used to linearise and discretise the convective term
        """
        ndim = self.simulation.ndim
        data = self.simulation.data
        
        # Update convective velocity field components
        for d in range(ndim):
            uic = data['u_conv%d' % d]
            uip =  data['up%d' % d]
            uipp = data['upp%d' % d]
            uic.vector()[:] = 2*uip.vector()[:] - uipp.vector()[:]
        
        # Update the convection blending factors
        for cs in self.convection_schemes:
            cs.update(t, dt, data['u_conv'])
    
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            a, L = self.eqs_mom_pred[d]
            
            # Solve the advection equation
            dolfin.solve(a == L, us, bcs)
    
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
        """
        p_hat = self.simulation.data['p_hat']
        bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        a, L = self.eq_pressure
        
        # Assemble matrix and vector
        A = dolfin.assemble(a)
        b = dolfin.assemble(L)
        
        # Create vector that spans the null space
        null_vec = dolfin.Vector(p_hat.vector())
        null_vec[:] = 1
        null_vec *= 1/null_vec.norm("l2")
        
        # Create null space basis object and attach to Krylov solver
        solver = dolfin.KrylovSolver(A, "gmres")
        null_space = dolfin.VectorSpaceBasis([null_vec])
        solver.set_nullspace(null_space)
        
        # Apply boundary conditions
        for bc in bcs:
            bc.apply(A, b)
        
        # Orthogonalize b with respect to the null space
        null_space.orthogonalize(b)
        
        solver.solve(p_hat.vector(), b) 
    
    def velocity_update(self):
        """
        Update the velocity estimates with the updated pressure
        field from the pressure correction equation
        """
        rho = self.simulation.data['rho']
        p_hat = self.simulation.data['p_hat']
        Vu = self.simulation.data['Vu']
        
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            u_new = self.simulation.data['u%d' % d]
            bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            
            # Update the velocity
            f = us - 2*self.dt/(3*rho) * p_hat.dx(d)
            un = dolfin.project(f, Vu, bcs)
            u_new.assign(un)
    
    def pressure_update(self):
        """
        Update the pressure at the end of an inner iteration
        """
        p_hat = self.simulation.data['p_hat']
        p = self.simulation.data['p']
        p.vector()[:] += p_hat.vector()[:] 
    
    def velocity_update_final(self):
        """
        Update the velocities at the end of the time step
        """
        for d in range(self.simulation.ndim):
            u_new = self.simulation.data['u%d' % d]
            us = self.simulation.data['u_star%d' % d]
            up = self.simulation.data['up%d' % d]
            upp = self.simulation.data['upp%d' % d]
            upp.vector()[:] = up.vector()[:]
            up.vector()[:] = u_new.vector()[:]
            us.vector()[:] = u_new.vector()[:]
    
    def run(self):
        """
        Run the simulation
        """
        dt = self.simulation.input['time']['dt']
        tmax = self.simulation.input['time']['tmax']
        num_inner_iter = self.simulation.input['solver'].get('num_inner_iter', 1)
        assert dt > 0
        
        t = 0
        it = 0
        while t+dt <= tmax + 1e-8:
            it += 1
            t += dt
            self.simulation.new_timestep(it, t, dt)
            self.dt.assign(dt)
            
            self.update_convection(t, dt)
            
            for _ in xrange(num_inner_iter):
                self.momentum_prediction(t, dt)
                self.pressure_correction()
                self.velocity_update()
                self.pressure_update()
            
            self.velocity_update_final()
            
            self.simulation.end_timestep()
