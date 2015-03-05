from __future__ import division
import dolfin
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error, timeit, velocity_error_norm
from . import Solver, register_solver
from .ipcs_equations import define_advection_problem, define_poisson_problem, define_penalty

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
        
        # Physical properties
        rho = sim.data['rho']
        nu = sim.data['nu']
        g = sim.data['g']
        self.dt = dolfin.Constant(1.0)
        
        # Get convection schemes for the velocity
        conv_schemes = []
        for d in range(sim.ndim):
            conv_scheme_name = sim.input.get_value('convection/u/convection_scheme', required_type='string')
            conv_scheme = get_convection_scheme(conv_scheme_name)(simulation, 'u%d' % d)
            conv_schemes.append(conv_scheme)
        self.convection_schemes = conv_schemes
        
        # Coefficients for u, up and upp 
        self.time_coeffs = dolfin.Constant([1, -1, 0]) # First time step
        self.relaxation_value = 1.0
        self.relaxation = dolfin.Constant(self.relaxation_value)
        
        # Define the momentum prediction equations
        Pu = Vu.ufl_element().degree()
        k = k_max = k_min = nu
        penalty = define_penalty(mesh, Pu, k_max, k_min)
        self.eqs_mom_pred = []
        for d in range(sim.ndim):
            trial = dolfin.TrialFunction(Vu)
            test = dolfin.TestFunction(Vu)
            beta = conv_schemes[d].blending_function
            f = -1/rho*p.dx(d) + g[d]
            dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            neumann_bcs = self.simulation.data['neumann_bcs'].get('u%d' % d, [])
            a1, L1 = define_advection_problem(trial, test, upvec[d], uppvec[d],
                                              u_conv, n, beta, self.time_coeffs, self.dt,
                                              dirichlet_bcs)
            a2, L2 = define_poisson_problem(trial, test, k, f, n, penalty, dirichlet_bcs, neumann_bcs)
            eq = a1+a2, L1+L2
            self.eqs_mom_pred.append(eq)
        
        # Define the pressure correction equation
        trial = dolfin.TrialFunction(Vp)
        test = dolfin.TestFunction(Vp)
        Pp = Vu.ufl_element().degree()
        k = k_max = k_min = 1/rho
        penalty = define_penalty(mesh, Pp, k_max, k_min)
        f = -self.time_coeffs[0]/self.dt * dolfin.nabla_div(u_star)
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        neumann_bcs = self.simulation.data['neumann_bcs'].get('p', [])
        self.eq_pressure = define_poisson_problem(trial, test, k, f, n, penalty, dirichlet_bcs, neumann_bcs)
        
        # For error norms in the convergence estimates
        elem = uvec[0].element()
        self.Vu_highp = dolfin.FunctionSpace(mesh, elem.family(), elem.degree() + 3)
    
    @timeit
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
            uic.vector()[:] = 5/3*uip.vector()[:] - 2/3*uipp.vector()[:]
            
            # Apply the Dirichlet boundary conditions
            if False:
                dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
                for bc in dirichlet_bcs:
                    bc.apply(uic.vector())
        
        # Update the convection blending factors
        for cs in self.convection_schemes:
            cs.update(t, dt, data['u'])
    
    @timeit
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            a, L = self.eqs_mom_pred[d]
            
            # Solve the advection equation
            family = us.element().family()
            if family == 'Lagrange':
                dolfin.solve(a == L, us, dirichlet_bcs)
            elif family == 'Discontinuous Lagrange':
                dolfin.solve(a == L, us)
                
                #for bc in dirichlet_bcs:
                #    bc.apply(us.vector())
    
    @timeit
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
        """
        p_hat = self.simulation.data['p_hat']
        dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('p', [])
        a, L = self.eq_pressure
        
        # Assemble matrix and vector
        A = dolfin.assemble(a)
        b = dolfin.assemble(L)
        
        # Apply strong boundary conditions
        family = p_hat.element().family()
        if family == 'Lagrange':
            for dbc in dirichlet_bcs:
                dbc.apply(A, b)
        
        # Get an appropriate solver
        solver = dolfin.KrylovSolver(A, "gmres")
        
        if not dirichlet_bcs:
            # Create vector that spans the null space
            null_vec = dolfin.Vector(p_hat.vector())
            null_vec[:] = 1
            null_vec *= 1/null_vec.norm("l2")
            
            # Create null space basis object and attach to Krylov solver
            null_space = dolfin.VectorSpaceBasis([null_vec])
            solver.set_nullspace(null_space)
            
            # Orthogonalize b with respect to the null space
            null_space.orthogonalize(b)
        
        solver.solve(p_hat.vector(), b) 
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity estimates with the updated pressure
        field from the pressure correction equation
        """
        rho = self.simulation.data['rho']
        p_hat = self.simulation.data['p_hat']
        Vu = self.simulation.data['Vu']
        c1 = self.time_coeffs[0]
        
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            u_new = self.simulation.data['u%d' % d]
            
            # Update the velocity
            f = us - self.relaxation*self.dt/(c1*rho) * p_hat.dx(d)
            un = dolfin.project(f, Vu)
            u_new.assign(un)
            
            # Apply the Dirichlet boundary conditions
            family = u_new.element().family()
            if family == 'Lagrange':
                # Apply strong boundary conditions
                dirichlet_bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
                for bc in dirichlet_bcs:
                    bc.apply(u_new.vector())
    
    @timeit
    def pressure_update(self):
        """
        Update the pressure at the end of an inner iteration
        """
        p_hat = self.simulation.data['p_hat']
        p = self.simulation.data['p']
        p.vector()[:] += self.relaxation_value*p_hat.vector()[:] 
    
    @timeit
    def velocity_update_final(self):
        """
        Update the velocities at the end of the time step
        """
        for d in range(self.simulation.ndim):
            u_new = self.simulation.data['u%d' % d]
            up = self.simulation.data['up%d' % d]
            upp = self.simulation.data['upp%d' % d]
            upp.vector()[:] = up.vector()[:]
            up.vector()[:] = u_new.vector()[:]
    
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation        
        t = 0
        it = 0
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            num_inner_iter = sim.input.get_value('solver/num_inner_iter', 1, required_type='int')
            allowable_error_inner = sim.input.get_value('solver/allowable_error_inner', 1e100, required_type='float')
            
            # Check if the simulation is done
            if t+dt > tmax + 1e-8:
                break
            
            # Advance one time step
            it += 1
            t += dt
            self.simulation.new_timestep(it, t, dt)
            self.dt.assign(dt)
            
            # Extrapolate the convecting velocity to the new time step
            self.update_convection(t, dt)
            
            for iit in xrange(num_inner_iter):
                self.momentum_prediction(t, dt)
                self.pressure_correction()
                self.velocity_update()
                self.pressure_update()
                
                # Convergence estimates
                L2s = velocity_error_norm(sim.data['u'], sim.data['u_star'], self.Vu_highp, 'L2')
                L2c = velocity_error_norm(sim.data['u'], sim.data['u_conv'], self.Vu_highp, 'L2')
                Linfs = velocity_error_norm(sim.data['u'], sim.data['u_star'], self.Vu_highp, 'Linf')
                Linfc = velocity_error_norm(sim.data['u'], sim.data['u_conv'], self.Vu_highp, 'Linf')
                sim.log.info('  Inner iteration %3d - Linf* %10.3e - Linfc %10.3e - L2* %10.3e - L2c %10.3e' % (iit+1, Linfs, Linfc, L2s, L2c))
                
                if Linfs < allowable_error_inner:
                    break
            
            self.velocity_update_final()
            
            # Change time coefficient to second order
            self.time_coeffs.assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Postprocess this time step
            sim.end_timestep()
        
        # We are done
        sim.end_simulation(success=True)
