import dolfin
from dolfin import dot, nabla_grad, avg, jump, dx, ds, dS
from ocellaris.convection import get_convection_scheme
from ocellaris.utils import report_error
from . import Solver, register_solver

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
        
        # Create pressure function
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['p_star'] = ps = dolfin.Function(Vp)
        
        # Mesh parameters
        n = dolfin.FacetNormal(mesh)
        h = dolfin.CellSize(mesh)
        
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
        penalty = 1.0/avg(h)
        self.eqs_mom_pred = []
        for d in range(sim.ndim):
            beta = conv_schemes[d].blending_function
            f = -1/rho*ps.dx(d) + g[d]
            a1, L1 = define_advection_problem(Vu, upvec[d], uppvec[d], 
                                              u_conv, n, beta, self.dt)
            a2, L2 = define_poisson_problem(Vu, nu, f, n, penalty, [])
            eq = a1+a2, L1+L2
            self.eqs_mom_pred.append(eq)
        
        # Define the pressure correction equation
        penalty = 1.0/avg(h)
        f = 1/rho
        neumann_bcs = self.simulation.data['neumann_bcs'].get('p', [])
        self.eq_pressure = define_poisson_problem(Vp, 1, f, n, penalty, neumann_bcs)
    
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            a, L = self.eqs_mom_pred[d]
            
            # Solve the advection equation
            self.dt.assign(dt)
            dolfin.solve(a == L, us, bcs)
    
    def pressure_correction(self):
        """
        Solve the pressure correction equation
        """
        pass
    
    def velocity_update(self):
        pass
    
    def run(self):
        """
        Run the simulation
        """
        ndim = self.simulation.ndim
        dt = self.simulation.input['time']['dt']
        tmax = self.simulation.input['time']['tmax']
        num_inner_iter = self.simulation.input['solver'].get('num_inner_iter', 1)
        assert dt > 0
        
        data = self.simulation.data
        t = 0
        it = 0
        while t+dt <= tmax + 1e-8:
            it += 1
            t += dt
            self.simulation.new_timestep(it, t, dt)
            self.dt.assign(dt)
            
            # Update convective velocity field components
            for d in range(ndim):
                uic = data['u_conv%d' % d]
                uip =  data['up%d' % d]
                uipp = data['upp%d' % d]
                uic.vector()[:] = 2*uip.vector()[:] - uipp.vector()[:] 
            
            # Update the convection blending factors
            for cs in self.convection_schemes:
                cs.update(t, dt, data['u_conv'])
            
            for _ in xrange(num_inner_iter):
                self.momentum_prediction(t, dt)
                self.pressure_correction()
                self.velocity_update()
            
            self.simulation.end_timestep()

def define_advection_problem(V, up, upp, u_conv, n, beta, dt):
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    
    # Upstream and downstream normal velocities
    flux_nU = u*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
    flux_nD = u*(dot(u_conv, n) - abs(dot(u_conv, n)))/2

    # Define the blended flux
    # The blending factor beta is not DG, so beta('+') == beta('-')
    b = beta('+')
    flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
    
    # Equation to solve
    eq = (3*u - 4*up + upp)/(2*dt)*v*dx \
         - u*dot(u_conv, nabla_grad(v))*dx \
         + flux*jump(v)*dS \
         + u*dot(u_conv, n)*v*ds
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)
    
    return a, L

def define_poisson_problem(V, k, f, n, penalty, neumann_bcs):
    """
    Define the Poisson problem for u in f.space V
    
        div(k*grad(u)) = f
    
    Returns the bilinear and linear forms
    """
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    
    # FIXME: introduce k in the equation!
    
    a = dot(nabla_grad(v), nabla_grad(u))*dx \
        - dot(avg(nabla_grad(v)), jump(u, n))*dS \
        - dot(jump(v, n), avg(nabla_grad(u)))*dS \
        + penalty*dot(jump(v, n), jump(u, n))*dS \
        - dot(nabla_grad(v), u*n)*ds \
        - dot(v*n, nabla_grad(u))*ds \
        + penalty*v*u*ds
    L = v*f*dx #- u0*dot(grad(v), n)*ds + (gamma/h)*u0*v*ds(1)
    
    # Add Neumann boundary conditions
    for nbc in neumann_bcs:
        L += v*nbc.value*nbc.ds 

    return a, L
