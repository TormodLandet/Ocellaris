# encoding: utf8
from __future__ import division
import dolfin
from dolfin import dx, div, grad, dot
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error, timeit, linear_solver_from_input
from . import Solver, register_solver, BDF, CRANK_NICOLSON


# Default values, can be changed in the input file
SOLVER = 'petsc'
LU_PARAMETERS = {}

# Implemented timestepping methods
TIMESTEPPING_METHODS = (BDF,)


@register_solver('Coupled')
class SolverIPCS(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the pressure-velocity splitting
        scheme IPCS (Incremental Pressure Correction Scheme)
        """
        self.simulation = sim = simulation
        self.create_functions()
        self.read_input()
        
        # First time step timestepping coefficients
        sim.data['time_coeffs'] = dolfin.Constant([1, -1, 0])
        self.is_first_timestep = True
        
        # Solver control parameters
        sim.data['dt'] = dolfin.Constant(simulation.dt)
        self.is_single_phase = isinstance(sim.data['rho'], dolfin.Constant)
        
        # Get the BCs for the coupled function space
        self.dirichlet_bcs = self.coupled_boundary_conditions()
        
        # Create equation
        self.define_coupled_equation()
        
        # Store number of iterations
        self.niters = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Solver for the coupled system
        self.coupled_solver = linear_solver_from_input(sim, 'solver/coupled', 'lu', 
                                                       None, SOLVER, LU_PARAMETERS)
        
        if isinstance(self.coupled_solver, dolfin.KrylovSolver):
            sim.log.warning('WARNING: Using a Krylov solver for the coupled NS equations is not a good idea')
        else:
            self.coupled_solver.parameters['same_nonzero_pattern'] = True
        
        # Coefficients for u, up and upp
        self.timestepping_method = sim.input.get_value('solver/timestepping_method', BDF, 'string')
        if not self.timestepping_method in TIMESTEPPING_METHODS:
            available_methods = '\n'.join(' - %s' % m for m in TIMESTEPPING_METHODS)
            report_error('Unknown timestepping method', 
                         'Timestepping method %s not recognised, please use one of:\n%s' %
                         (self.timestepping_method, available_methods))
    
    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation
        
        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        
        if sim.ndim == 2:
            self.Vcoupled = Vu * Vu * Vp
            self.coupled_func = dolfin.Function(self.Vcoupled)
            
            self.subspaces = [self.Vcoupled.sub(0).sub(0), self.Vcoupled.sub(0).sub(1),
                              self.Vcoupled.sub(1)]
            #self.subspace_vectors = [self.coupled_func.sub(0)[0].vector(),
            #                         self.coupled_func.sub(0)[1].vector(),
            #                         self.coupled_func.sub(1).vector()]
            self.subspace_names = 'u0 u1 p'.split()
        
        else:
            self.Vcoupled = Vu * Vu * Vu * Vp
            self.coupled_func = dolfin.Function(self.Vcoupled)
            
            self.subspaces = [self.Vcoupled.sub(0).sub(0), self.Vcoupled.sub(0).sub(1),
                              self.Vcoupled.sub(0).sub(2), self.Vcoupled.sub(1)]
            #self.subspace_vectors = [self.coupled_func.sub(0)[0].vector(),
            #                         self.coupled_func.sub(0)[1].vector(),
            #                         self.coupled_func.sub(0)[2].vector(),
            #                         self.coupled_func.sub(1).vector()]
            self.subspace_names = 'u0 u1 u2 p'.split()
        
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
        sim.data['coupled'] = dolfin.Function
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        
    def define_coupled_equation(self):
        """
        Setup the coupled Navier-Stokes equation
        """
        u, p = dolfin.TrialFunctions(self.Vcoupled)
        v, q = dolfin.TestFunctions(self.Vcoupled)
        
        c1, c2, c3 = self.simulation.data['time_coeffs']
        dt = self.simulation.data['dt']
        u_conv = self.simulation.data['u_conv']
        u_expl = self.simulation.data['up']
        g = self.simulation.data['g']
        rho = self.simulation.data['rho']
        nu = self.simulation.data['nu']
        mu = rho*nu
        
        div, grad = dolfin.nabla_div, dolfin.nabla_grad
        
        # Divergence free criterion
        # ∇⋅u = 0
        eq = -div(u)*q*dx
        
        # Momentum equations
        for d in range(self.simulation.ndim):
            up = self.simulation.data['up%d' % d]
            upp = self.simulation.data['upp%d' % d]
            
            # Time derivative
            # ∂u/∂t
            eq += rho*(c1*u[d] + c2*up + c3*upp)/dt*v[d]*dx
            
            # Convection
            # ∇⋅(ρ u ⊗ u_conv)
            eq += div(rho*u[d]*u_conv)*v[d]*dx
            #eq += rho*v[d]*dot(u_conv, grad(u[d]))*dx
            
            # Diffusion
            # -∇⋅μ[(∇u) + (∇u_expl)^T]
            eq += mu*dot(grad(u[d]), grad(v[d]))*dx
            eq += mu*dot(u_expl.dx(d), grad(v[d]))*dx
            
            # Pressure
            # ∇p
            eq -= v[d].dx(d)*p*dx
            
            # Body force (gravity)
            # ρ g
            eq -= rho*g[d]*v[d]*dx
        
        self.eq_coupled = dolfin.system(eq)
    
    def coupled_boundary_conditions(self):
        """
        Convert boundary conditions from segregated to coupled function spaces
        """
        coupled_dirichlet_bcs = []
        for i, name in enumerate(self.subspace_names):
            V = self.subspaces[i]
            bcs = self.simulation.data['dirichlet_bcs'].get(name, [])
            for bc in bcs:
                bc_new = bc.copy_and_change_function_space(V)
                coupled_dirichlet_bcs.append(bc_new)
        
        return coupled_dirichlet_bcs
    
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
                uic.vector()[:] = uip.vector()[:]
            else:
                # Backwards difference formulation - standard linear extrapolation 
                uic.vector()[:] = 2.0*uip.vector()[:] - 1.0*uipp.vector()[:]
        
        self.is_first_timestep = False
    
    @timeit
    def solve_coupled(self):
        """
        Solve the coupled equations
        """
        # Assemble the equation system
        a, L = self.eq_coupled
        A, b = dolfin.assemble_system(a, L)
        
        # Apply strong boundary conditions
        for dbc in self.dirichlet_bcs:
            dbc.apply(A, b)
        
        # Solve the equation system
        self.coupled_solver.solve(A, self.coupled_func.vector(), b)
        
        # Copy the data from the coupled solution vector and to the individual components
        u, p = self.coupled_func.split(True)
        funcs = list(u.split(True)) + [p]
        for name, func in zip(self.subspace_names, funcs):
            self.simulation.data[name].vector()[:] = func.vector()
    
    @timeit
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
            this_maxabs = abs(sim.data['upp%d' % d].vector().array()).max()
            maxabs = max(maxabs, this_maxabs)
        has_upp_start_values = maxabs > 0 
        
        # Previous-previous values are provided so we can start up with second order time stepping 
        if has_upp_start_values:
            sim.log.info('Initial values for upp are found and used')
            self.is_first_timestep = False
            if self.timestepping_method == BDF:
                self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
        
        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')
            
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
            
            # Solve for the new time step
            self.solve_coupled()
            
            # Move u -> up, up -> upp and prepare for the next time step
            for d in range(self.simulation.ndim):
                u_new = self.simulation.data['u%d' % d]
                up = self.simulation.data['up%d' % d]
                upp = self.simulation.data['upp%d' % d]
                upp.assign(up)
                up.assign(u_new)
            
            # Change time coefficient to second order
            if self.timestepping_method == BDF:
                self.simulation.data['time_coeffs'].assign(dolfin.Constant([3/2, -2, 1/2]))
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
