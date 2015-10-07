# encoding: utf8
from __future__ import division
import dolfin
from dolfin import Constant
from ocellaris.utils import report_error, timeit, linear_solver_from_input
from . import Solver, register_solver, BDF, BDM
from .coupled_equations import CoupledEquations, LAGRANGE_MULTIPLICATOR
from .dg_helpers import VelocityBDMProjection


# Default values, can be changed in the input file
LU_SOLVER = 'petsc'
LU_PARAMETERS = {}

# Implemented timestepping methods
TIMESTEPPING_METHODS = (BDF,)


@register_solver('Coupled')
class SolverCoupled(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the pressure-velocity splitting
        scheme IPCS (Incremental Pressure Correction Scheme)
        """
        self.simulation = sim = simulation
        self.read_input()
        self.create_functions()
        
        # First time step timestepping coefficients
        self.set_timestepping_coefficients([1, -1, 0])
        
        # Solver control parameters
        sim.data['dt'] = Constant(simulation.dt)
        self.is_single_phase = isinstance(sim.data['rho'], Constant)
        
        # Get the BCs for the coupled function space
        self.dirichlet_bcs = self.coupled_boundary_conditions()
        
        # Create equation
        self.eqs = CoupledEquations(simulation, self.timestepping_method)
        
        # Velocity post_processing
        self.velocity_postprocessor = None
        if self.velocity_postprocessing_method == BDM:
            self.velocity_postprocessor = VelocityBDMProjection(sim.data['u'])
        
        
        # Store number of iterations
        self.niters = None
    
    def read_input(self):
        """
        Read the simulation input
        """
        sim = self.simulation
        
        # Solver for the coupled system
        self.coupled_solver = linear_solver_from_input(sim, 'solver/coupled', 'lu', 
                                                       None, LU_SOLVER, LU_PARAMETERS)
        
        # Give warning if using iterative solver
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
        
        # Lagrange multiplicator or remove null space via PETSc or just normalize after solving
        self.remove_null_space = True
        self.pressure_null_space = None
        self.use_lagrange_multiplicator = sim.input.get_value('solver/use_lagrange_multiplicator',
                                                              LAGRANGE_MULTIPLICATOR, 'bool')
        if self.use_lagrange_multiplicator or self.simulation.data['dirichlet_bcs'].get('p', []):
            self.remove_null_space = False
        
        # Check if the solver supports removing null spaces
        self.normalize_pressure = False
        does_not_support_null_space = ('mumps', )
        if self.remove_null_space and self.coupled_solver.created_with_lu_method in does_not_support_null_space:    
            self.normalize_pressure = True
            self.remove_null_space = False
            
        # Representation of velocity
        Vu_family = sim.data['Vu'].ufl_element().family()
        self.vel_is_discontinuous = (Vu_family == 'Discontinuous Lagrange')
        
        # Local DG velocity postprocessing 
        default_postprocessing = BDM if self.vel_is_discontinuous else None
        self.velocity_postprocessing_method = sim.input.get_value('solver/velocity_postprocessing',
                                                                  default_postprocessing, 'string')
    
    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation
        
        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        cd = sim.data['constrained_domain']
            
        # Create coupled mixed function space and mixed function to hold results
        func_spaces = [Vu] * sim.ndim + [Vp]
        self.subspace_names = ['u%d' % d for d in range(sim.ndim)] + ['p']
        
        if self.use_lagrange_multiplicator:
            Vl = dolfin.FunctionSpace(sim.data['mesh'], "R", 0, constrained_domain=cd)
            sim.data['l'] = dolfin.Function(Vl)
            func_spaces.append(Vl)
            self.subspace_names.append('l')
        
        Vcoupled = dolfin.MixedFunctionSpace(func_spaces)
        sim.data['Vcoupled'] = Vcoupled
        
        Nspace = len(func_spaces)
        self.subspaces = [Vcoupled.sub(i) for i in range(Nspace)]
        sim.data['coupled'] = self.coupled_func = dolfin.Function(Vcoupled)
        self.assigner = dolfin.FunctionAssigner(func_spaces, Vcoupled)
        
        # Create segregated functions on component and vector form
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
        sim.data['p'] = dolfin.Function(Vp)
    
    def coupled_boundary_conditions(self):
        """
        Convert boundary conditions from segregated to coupled function spaces
        """
        coupled_dirichlet_bcs = []
        for i, name in enumerate(self.subspace_names):
            if self.vel_is_discontinuous and name.startswith('u'):
                # Use weak BCs if the velocity is DG
                continue
            
            V = self.subspaces[i]
            bcs = self.simulation.data['dirichlet_bcs'].get(name, [])
            for bc in bcs:
                bc_new = bc.copy_and_change_function_space(V)
                coupled_dirichlet_bcs.append(bc_new)
        
        return coupled_dirichlet_bcs

    def set_timestepping_coefficients(self, coeffs):
        """
        Set the time stepping coefficients used for the temporal derivative
        """
        if not 'time_coeffs' in self.simulation.data:
            self.is_first_timestep = True 
            self.simulation.data['time_coeffs'] = Constant(coeffs)
            self.simulation.data['time_coeffs_py'] = coeffs
        else:
            self.simulation.data['time_coeffs'].assign(Constant(coeffs))
            self.simulation.data['time_coeffs_py'] = coeffs
    
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
    def postprocess_velocity(self):
        """
        Apply a post-processing operator to the given velocity field
        """
        if self.velocity_postprocessor:
            self.velocity_postprocessor.run()
    
    @timeit
    def solve_coupled(self):
        """
        Solve the coupled equations
        """
        # Assemble the equation system
        A = self.eqs.assemble_lhs()
        b = self.eqs.assemble_rhs()
        
        # Apply strong boundary conditions (this list is empty for DG)
        for dbc in self.dirichlet_bcs:
            dbc.apply(A, b)
        
        if self.remove_null_space:
            if self.pressure_null_space is None:
                # Create null space vector in Vp Space
                null_func = dolfin.Function(self.simulation.data['Vp'])
                null_vec = null_func.vector()
                null_vec[:] = 1
                null_vec *= 1/null_vec.norm("l2")
                
                # Convert null space vector to coupled space
                null_func2 = dolfin.Function(self.simulation.data['Vcoupled'])
                ndim = self.simulation.ndim
                fa = dolfin.FunctionAssigner(self.subspaces[ndim], self.simulation.data['Vp'])
                fa.assign(null_func2.sub(ndim), null_func)
                
                # Create the null space basis
                self.pressure_null_space = dolfin.VectorSpaceBasis([null_func2.vector()])
            
            # Make sure the null space is set on the matrix
            dolfin.as_backend_type(A).set_nullspace(self.pressure_null_space)
            
            # Orthogonalize b with respect to the null space
            self.pressure_null_space.orthogonalize(b)
        
        # Solve the equation system
        self.coupled_solver.solve(A, self.coupled_func.vector(), b)
        
        # Assign into the regular (split) functions from the coupled function
        funcs = [self.simulation.data[name] for name in self.subspace_names]
        self.assigner.assign(funcs, self.coupled_func)
        for func in funcs:
            func.vector().apply('insert') # dolfin bug #587
                    
        # Some solver cannot remove the null space, so we just normalize the pressure instead
        if self.normalize_pressure:
            p = self.simulation.data['p']
            dx2 = dolfin.dx(domain=p.function_space().mesh())
            vol = dolfin.assemble(dolfin.Constant(1)*dx2)
            # Perform correction twice due to round-of error. The first correction can be
            # i.e 1e14 while the next correction is around unity
            for _ in range(2):
                pavg = dolfin.assemble(p*dx2)
                p.vector()[:] -= pavg/vol   
    
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
                self.set_timestepping_coefficients([3/2, -2, 1/2])
        
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
            
            # Postprocess the solution velocity field
            self.postprocess_velocity()
            
            # Move u -> up, up -> upp and prepare for the next time step
            for d in range(self.simulation.ndim):
                u_new = self.simulation.data['u%d' % d]
                up = self.simulation.data['up%d' % d]
                upp = self.simulation.data['upp%d' % d]
                upp.assign(up)
                up.assign(u_new)
            
            # Change time coefficient to second order
            if self.timestepping_method == BDF:
                self.set_timestepping_coefficients([3/2, -2, 1/2])
            
            # Postprocess this time step
            sim.hooks.end_timestep()
        
        # We are done
        sim.hooks.simulation_ended(success=True)
