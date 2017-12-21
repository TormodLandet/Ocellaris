import dolfin as df
from dolfin import dot, sqrt, grad, jump, avg, dx, dS, Form
from ocellaris.utils import timeit


class SolutionProperties(object):
    def __init__(self, simulation):
        """
        Calculate Courant and Peclet numbers
        """
        self.simulation = simulation
        self.divergence_method = None
        self.active = False
        self.has_div_conv = False
        self._div = {}
    
    def setup(self):
        sim = self.simulation
        
        self.active = sim.input.get_value('output/solution_properties', True, 'bool')
        divergence_method = sim.input.get_value('output/divergence_method', 'div', 'string')
        plot_divergences = sim.input.get_value('output/plot_divergences', False, 'bool')
        if not self.active:
            sim.log.info('SolutionProperties not active')
            return
        sim.log.info('SolutionProperties active with div method %r' % divergence_method)
        
        self.mesh = sim.data['mesh']
        u = sim.data['u']
        dt = sim.data['dt']
        rho = sim.data['rho']
        nu = sim.data['nu']
        g = sim.data['g']
        x0 = df.Constant([0]*self.mesh.topology().dim())
    
        self._setup_courant(u, dt)
        self._setup_peclet(u, nu)
        self._setup_divergence(u, divergence_method)
        self._setup_energy(rho, u, g, x0)
        self._setup_mass(rho)
        
        if 'up_conv0' in sim.data:
            u_conv = sim.data['u_conv']
            self._setup_divergence(u_conv, divergence_method, 'u_conv')
            self.has_div_conv = True
        
        if plot_divergences:
            for name in ('u', 'u_conv'):
                if name in self._div:
                    sim.io.add_extra_output_function(self._div[name]['div_dS'])
                    sim.io.add_extra_output_function(self._div[name]['div_dx'])
    
    def _setup_courant(self, vel, dt):
        """
        Co = a*dt/h where a = mag(vel)
        """
        V = df.FunctionSpace(self.mesh, 'DG', 0)
        h = self.simulation.data['h']
        
        u, v = df.TrialFunction(V), df.TestFunction(V)
        a = u*v*dx
        vmag = sqrt(dot(vel, vel))
        L = vmag*dt/h*v*dx
        
        # Pre-factorize matrices and store for usage in projection
        self._courant_solver = df.LocalSolver(a, L)
        self._courant_solver.factorize()
        self._courant = df.Function(V)
    
    def _setup_peclet(self, vel, nu):
        """
        Pe = a*h/(2*nu) where a = mag(vel)
        """
        V = df.FunctionSpace(self.mesh, 'DG', 0)
        h = self.simulation.data['h']
        u, v = df.TrialFunction(V), df.TestFunction(V)
        a = u*v*dx
        L = dot(vel, vel)**0.5*h/(2*nu)*v*dx
        
        # Pre-factorize matrices and store for usage in projection
        self._peclet_solver = df.LocalSolver(a, L)
        self._peclet_solver.factorize()
        self._peclet = df.Function(V)
    
    def _setup_divergence(self, vel, method, name='u'):
        """
        Calculate divergence and element to element velocity
        flux differences on the same edges
        """
        V = df.FunctionSpace(self.mesh, 'DG', 0)
        n = df.FacetNormal(self.mesh)
        u, v = df.TrialFunction(V), df.TestFunction(V)
        
        # The difference between the flux on the same facet between two different cells
        a1 = u*v*dx
        w = dot(vel('+') - vel('-'), n('+'))
        if method == 'div0':
            L1 = w*avg(v)*dS
        else:
            L1 = abs(w)*avg(v)*dS
        
        # The divergence internally in the cell
        a2 = u*v*dx
        if method in ('div', 'div0'):
            L2 = abs(df.div(vel))*v*dx
        elif method == 'gradq_avg':
            L2 = dot(avg(vel), n('+'))*jump(v)*dS - dot(vel, grad(v))*dx
        else:
            raise ValueError('Divergence type %r not supported' % method)
        
        # Store for usage in projection
        storage = self._div[name] = {}
        storage['dS_solver'] = df.LocalSolver(a1, L1)
        storage['dx_solver'] = df.LocalSolver(a2, L2)
        storage['div_dS'] = df.Function(V)
        storage['div_dx'] = df.Function(V)
        
        # Pre-factorize matrices
        storage['dS_solver'].factorize()
        storage['dx_solver'].factorize()
        storage['div_dS'].rename('Divergence_%s_dS' % name, 'Divergence_%s_dS' % name)
        storage['div_dx'].rename('Divergence_%s_dx' % name, 'Divergence_%s_dx' % name)
    
    def _setup_energy(self, rho, vel, gvec, x0):
        """
        Calculate kinetic and potential energy
        """
        x = df.SpatialCoordinate(self.mesh)
        self._form_E_k = Form(1/2*rho*dot(vel, vel)*dx(domain=self.mesh))
        self._form_E_p = Form(rho*dot(-gvec, x - x0)*dx(domain=self.mesh))
    
    def _setup_mass(self, rho):
        """
        Calculate mass
        """
        self._form_mass = Form(rho*dx(domain=self.mesh))
    
    @timeit
    def courant_number(self):
        """
        Calculate the Courant numbers in each cell
        """
        self._courant_solver.solve_local_rhs(self._courant)
        return self._courant
    
    @timeit
    def peclet_number(self):
        """
        Calculate the Peclet numbers in each cell
        """
        self._peclet_solver.solve_local_rhs(self._peclet)
        return self._peclet
    
    @timeit
    def divergences(self, name='u'):
        """
        Calculate the difference between the flux on the same
        facet between two different cells and the divergence
        inside each cell.
        
        Returns the sum of facet errors for each cell and the
        divergence error in each cell as DG0 functions 
        """
        storage = self._div[name]
        storage['dS_solver'].solve_global_rhs(storage['div_dS'])
        storage['dx_solver'].solve_global_rhs(storage['div_dx'])
        return storage['div_dS'], storage['div_dx']
    
    @timeit
    def total_energy(self):
        """
        Calculate the total energy in the field
        """
        E_k = df.assemble(self._form_E_k)
        E_p = df.assemble(self._form_E_p)
        return E_k, E_p
    
    @timeit
    def total_mass(self):
        """
        Calculate the total mass
        """
        mass = df.assemble(self._form_mass)
        return mass
