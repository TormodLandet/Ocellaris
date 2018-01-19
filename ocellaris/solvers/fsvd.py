import numpy
import dolfin
from ocellaris.utils import verify_key, timeit, linear_solver_from_input, ocellaris_error
from . import Solver, register_solver
from ..solver_parts import VelocityBDMProjection, SlopeLimiter
from .fsvd_equations import EQUATION_SUBTYPES


# Default values, can be changed in the input file
SOLVER_U = 'gmres'
PRECONDITIONER_U = 'additive_schwarz'
SOLVER_P = 'gmres'
PRECONDITIONER_P = 'hypre_amg'
KRYLOV_PARAMETERS = {'nonzero_initial_guess': True,
                     'relative_tolerance': 1e-10,
                     'absolute_tolerance': 1e-15}
EQUATION_SUBTYPE = 'Default'
USE_GRAD_P_FORM = False
FIX_PRESSURE_DOF = False


@register_solver('FSVD')
class SolverFSVD(Solver):
    def __init__(self, simulation):
        """
        A fractional step variable density incompressible Navier-Stokes solver
        """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        self.is_first_timestep = True
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        
        # Get equations
        MomentumPredictionEquation, PressureCorrectionEquation, \
            VelocityUpdateEquation = EQUATION_SUBTYPES[self.equation_subtype]
        
        # Equations to solve
        self.eqs_mom_pred = []
        for d in range(sim.ndim):
            eq = MomentumPredictionEquation(simulation, d, self.use_grad_p_form)
            self.eqs_mom_pred.append(eq)
        self.eq_pressure = PressureCorrectionEquation(simulation)
        self.eqs_vel_upd = []
        for d in range(sim.ndim):
            eq = VelocityUpdateEquation(simulation, d)
            self.eqs_vel_upd.append(eq)
        self.velocity_postprocessor = VelocityBDMProjection(sim, sim.data['u'])
        
        # Slope limiter for the momenum equation velocity components
        self.slope_limiters = [SlopeLimiter(sim, 'u', sim.data['u%d' % d], 'u%d' % d, old_value=sim.data['up%d' % d])
                               for d in range(sim.ndim)]
        
        # Pre assembled matrices
        self.Au = [None]*sim.ndim
        self.Ap = None
        self.pressure_null_space = None
        self.pressure_row_to_fix = numpy.array([0], dtype=numpy.intc)
        
        # Store the number of iterations used in the linear solvers
        self.niters_u = [None] * sim.ndim
        self.niters_p = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Check the representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        Vp_family = sim.data['Vp'].ufl_element().family()
        verify_key('finite element family', Vu_family, ('Discontinuous Lagrange',), 'FSCD velocity solver')
        verify_key('finite element family', Vp_family, ('Discontinuous Lagrange',), 'FSCD pressure solver')
        
        # Create linear solvers
        self.velocity_solver = linear_solver_from_input(self.simulation, 'solver/u', SOLVER_U,
                                                        PRECONDITIONER_U, None, KRYLOV_PARAMETERS)
        self.pressure_solver = linear_solver_from_input(self.simulation, 'solver/p', SOLVER_P,
                                                        PRECONDITIONER_P, None, KRYLOV_PARAMETERS)
        self.u_upd_solver = None
        
        # Get the class to be used for the equation system assembly
        self.equation_subtype = sim.input.get_value('solver/equation_subtype', EQUATION_SUBTYPE, 'string')
        verify_key('equation sub-type', self.equation_subtype, EQUATION_SUBTYPES, 'fsvd solver')
        
        # No need for special treatment if the pressure is set via Dirichlet conditions somewhere
        self.remove_null_space = len(sim.data['dirichlet_bcs'].get('p', [])) == 0
        
        # Control the form of the governing equations 
        self.use_grad_p_form = sim.input.get_value('solver/use_grad_p_form', USE_GRAD_P_FORM, 'bool')
        self.fix_pressure_dof = sim.input.get_value('solver/fix_pressure_dof', FIX_PRESSURE_DOF, 'bool')
    
    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation
        
        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        
        # Create velocity functions. Keep both component and vector forms
        u_list, up_list, upp_list, u_conv = [], [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            sim.data['u_conv%d' % d] = uc = dolfin.Function(Vu)
            u_list.append(u)
            up_list.append(up)
            upp_list.append(upp)
            u_conv.append(uc)
        sim.data['u'] = dolfin.as_vector(u_list)
        sim.data['up'] = dolfin.as_vector(up_list)
        sim.data['upp'] = dolfin.as_vector(upp_list)
        sim.data['u_conv'] = dolfin.as_vector(u_conv)
        self.u_tmp = dolfin.Function(Vu)
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['p_hat'] = dolfin.Function(Vp)
    
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
            
            if self.is_first_timestep:
                uic.assign(uip)
            else:
                uic.vector().zero()
                uic.vector().axpy(2.0, uip.vector())
                uic.vector().axpy(-1.0, uipp.vector())
        
        self.is_first_timestep = False
    
    @timeit
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        solver = self.velocity_solver
        
        err = 0.0
        for d in range(self.simulation.ndim):
            u = self.simulation.data['u%d' % d]
            self.u_tmp.assign(u)
            
            eq = self.eqs_mom_pred[d]
            
            if self.inner_iteration == 1:
                # Assemble the A matrix only the first inner iteration
                self.Au[d] = eq.assemble_lhs()
            
            A = self.Au[d]
            b = eq.assemble_rhs()
            
            self.simulation.hooks.matrix_ready('Au%d' % d, A, b)
            self.niters_u[d] = solver.solve(A, u.vector(), b)
            self.slope_limiters[d].run()
            
            self.u_tmp.vector().axpy(-1, u.vector())
            err += self.u_tmp.vector().norm('l2')
        return err
    
    @timeit
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        
        We handle the case where only Neumann conditions are given
        for the pressure by taking out the nullspace, a constant shift
        of the pressure, by providing the nullspace to the solver
        """
        p = self.simulation.data['p']
        
        # Assemble the A matrix only the first inner iteration
        if self.Ap is None:
            self.Ap = self.eq_pressure.assemble_lhs()
        
        # The equation system to solve
        A = self.Ap
        b = self.eq_pressure.assemble_rhs()
        
        if self.fix_pressure_dof:
            A.ident(self.pressure_row_to_fix)
        elif self.remove_null_space:  
            if self.pressure_null_space is None:
                null_vec = dolfin.Vector(p.vector())
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_vec])
                dolfin.as_backend_type(A).set_nullspace(self.pressure_null_space)
            self.pressure_null_space.orthogonalize(b)
    
        # Temporarily store the old pressure
        p_hat = self.simulation.data['p_hat']
        p_hat.vector().zero()
        p_hat.vector().axpy(-1, p.vector())
        
        # Solve for new pressure
        self.simulation.hooks.matrix_ready('Ap', A, b)
        self.niters_p = self.pressure_solver.solve(A, p.vector(), b)
        
        # Removing the null space of the matrix system is not strictly the same as removing
        # the null space of the equation, so we correct for this here 
        if self.remove_null_space:
            dx2 = dolfin.dx(domain=p.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1)*dx2)
            pavg = dolfin.assemble(p*dx2)/vol
            p.vector()[:] -= pavg
        
        # Calculate p_hat = p_new - p_old 
        p_hat.vector().axpy(1, p.vector())
        
        return p_hat.vector().norm('l2')
    
    @timeit
    def velocity_update(self):
        """
        Update the velocity predictions with the updated pressure
        field from the pressure correction equation
        """
        if self.u_upd_solver is None:
            self.u_upd_solver = dolfin.LocalSolver(self.eqs_vel_upd[0].form_lhs)
            self.u_upd_solver.factorize()
        
        Vu = self.simulation.data['Vu']
        for d in range(self.simulation.ndim):
            eq = self.eqs_vel_upd[d]
            b = eq.assemble_rhs()
            u = self.simulation.data['u%d' % d]
            self.u_upd_solver.solve_local(u.vector(), b, Vu.dofmap())
    
    @timeit
    def postprocess_velocity(self):
        """
        Apply a post-processing operator to the given velocity field
        """
        self.velocity_postprocessor.run()
    
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation        
        sim.hooks.simulation_started()
        t = sim.time
        it = sim.timestep
        
        # Check if there are non-zero values in the upp vectors
        maxabs = 0
        for d in range(sim.ndim):
            this_maxabs = abs(sim.data['upp%d' % d].vector().get_local()).max()
            maxabs = max(maxabs, this_maxabs)
        maxabs = dolfin.MPI.max(dolfin.MPI.comm_world, float(maxabs))
        has_upp_start_values = maxabs > 0
        
        # Previous-previous values are provided so we can start up with second order time stepping 
        if has_upp_start_values:
            sim.log.info('Initial values for upp are found and used')
            self.is_first_timestep = False
            self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
        
        # Give reasonable starting guesses for the solvers
        for d in range(sim.ndim):
            up = self.simulation.data['up%d' % d]
            u_new = self.simulation.data['u%d' % d]
            u_new.assign(up)
        
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            num_inner_iter = sim.input.get_value('solver/num_inner_iter', 100, 'int')
            allowable_div_inner = sim.input.get_value('solver/allowable_div_inner', 1e-6, 'float')
            
            # Check if the simulation is done
            if t+dt > tmax + 1e-6:
                break
            
            # Advance one time step
            it += 1
            t += dt
            self.simulation.data['dt'].assign(dt)
            self.simulation.hooks.new_timestep(it, t, dt)
            
            # Extrapolate the convecting velocity to the new time step
            self.update_convection(t, dt)
            
            # Run inner iterations
            self.inner_iteration = 1
            while self.inner_iteration <= num_inner_iter:
                err_u_star = self.momentum_prediction(t, dt)
                err_p = self.pressure_correction()
                
                div_dS_f, div_dx_f = sim.solution_properties.divergences()
                err_div = div_dS_f.vector().max() + div_dx_f.vector().max()
                
                # Information from solvers regarding number of iterations needed to solve linear system
                niters = ['%3d u%d' % (ni, d) for d, ni in enumerate(self.niters_u)]
                niters.append('%3d p' % self.niters_p)
                solver_info = ' - iters: %s' % ' '.join(niters)

                # Get max u_star
                umax = 0
                for d in range(sim.ndim):
                    thismax = abs(sim.data['u%d' % d].vector().get_local()).max()
                    umax = max(thismax, umax)
                umax = dolfin.MPI.max(dolfin.MPI.comm_world, float(umax))
                
                # Convergence estimates
                sim.log.info('  Inner iteration %3d - err u* %10.3e - err p %10.3e%s  ui*max %10.3e'
                             % (self.inner_iteration, err_u_star, err_p, solver_info,  umax)
                             + ' err div %10.3e' % err_div)
                
                if err_div < allowable_div_inner:
                    break
                elif self.inner_iteration > 3 and err_div > 1e4:
                    ocellaris_error('Iteration diverged', 'Inner iterations diverged')
                
                self.inner_iteration += 1
            
            self.velocity_update()
            self.postprocess_velocity()
            
            # Move u -> up, up -> upp and prepare for the next time step
            for d in range(sim.ndim):
                u_new = sim.data['u%d' % d]
                up = sim.data['up%d' % d]
                upp = sim.data['upp%d' % d]
                upp.assign(up)
                up.assign(u_new)
            
            # Change time coefficient to second order
            sim.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
