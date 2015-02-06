import dolfin
from dolfin import dot, nabla_grad, avg, jump, dx, ds, dS
from ocellaris.convection import get_convection_scheme
from . import Solver, register_solver
from scipy import ndimage

@register_solver('IPCS')
class SolverIPCS(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes solver based on the pressure-velocity splitting
        scheme IPCS (Incremental Pressure Correction Scheme)
        """
        self.simulation = sim = simulation
        mesh = sim.data['mesh']
        
        # Construct function spaces
        Pu = sim.input['solver'].get('polynomial_degree_velocity', 1)
        Pp = sim.input['solver'].get('polynomial_degree_pressure', 1)
        Vu = dolfin.FunctionSpace(mesh, 'Discontinuous Lagrange', Pu)
        Vp = dolfin.FunctionSpace(mesh, 'Discontinuous Lagrange', Pp)
        
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
            a1, L1 = define_advection_problem(Vu, u_conv, n, beta, self.dt)
            a2, L2 = define_poisson_problem(Vu, f, n, penalty)
            eq = a1+a2, L1+L2
            self.eqs_mom_pred.append(eq)
        
        # Define the pressure correction equation
        penalty = 1.0/avg(h)
        f = dolfin.Constant(2)
        self.eq_pressure = define_poisson_problem(Vp, f, n, penalty)
    
    def momentum_prediction(self, t, dt):
        """
        Solve the momentum prediction equation
        """
        for d in range(self.simulation.ndim):
            us = self.simulation.data['u_star%d' % d]
            bcs = self.simulation.data['dirichlet_bcs'].get('u%d' % d, [])
            a, L = self.eqs_mom_pred[d]
            bc = [] # FIXME
            
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

def define_advection_problem(V, u_conv, n, beta, dt):
    return 0, 0

def define_poisson_problem(V, f, n, penalty):
    """
    Define the Poisson problem
    
        div(grad(u)) = f
    
    Returns the bilinear and linear forms
    """
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    a = dot(nabla_grad(v), nabla_grad(u))*dx \
        - dot(avg(nabla_grad(v)), jump(u, n))*dS \
        - dot(jump(v, n), avg(nabla_grad(u)))*dS \
        + penalty*dot(jump(v, n), jump(u, n))*dS \
        - dot(nabla_grad(v), u*n)*ds \
        - dot(v*n, nabla_grad(u))*ds \
        + penalty*v*u*ds
    L = v*f*dx #- u0*dot(grad(v), n)*ds + (gamma/h)*u0*v*ds(1) + g*v*ds(2)

    return a, L
