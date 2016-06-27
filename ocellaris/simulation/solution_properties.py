# enconding: utf8
from __future__ import division
import dolfin as df
from dolfin import dot, sqrt, grad, jump, avg, dx, dS
from ocellaris.utils import timeit


class SolutionProperties(object):
    def __init__(self, simulation, divergence='div'):
        """
        Calculate Courant and Peclet numbers
        """
        self.simulation = simulation
        self.divergence_method = divergence
    
    def setup(self):
        sim = self.simulation
        self.mesh = sim.data['mesh']
        
        u = sim.data['u']
        dt = sim.data['dt']
        rho = sim.data['rho']
        nu = sim.data['nu']
        g = sim.data['g']
        x0 = df.Constant([0]*self.mesh.topology().dim())
    
        self._setup_courant(u, dt)
        self._setup_peclet(u, nu)
        self._setup_divergence(u, self.divergence_method)
        self._setup_energy(rho, u, g, x0)
        self._setup_mass(rho)
    
    def _setup_courant(self, vel, dt):
        """
        Co = a*dt/h where a = mag(vel)
        """
        V = df.FunctionSpace(self.mesh, 'DG', 0)
        h = df.CellSize(self.mesh)
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
        h = df.CellSize(self.mesh)
        df_nu = df.Constant(nu)
        u, v = df.TrialFunction(V), df.TestFunction(V)
        a = u*v*dx
        L = dot(vel, vel)**0.5*h/(2*df_nu)*v*dx
        
        # Pre-factorize matrices and store for usage in projection
        self._peclet_solver = df.LocalSolver(a, L)
        self._peclet_solver.factorize()
        self._peclet = df.Function(V)
    
    def _setup_divergence(self, vel, method):
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
        L1 = abs(w)*avg(v)*dS
        
        # The divergence internally in the cell
        a2 = u*v*dx
        if method == 'div':
            L2 = abs(df.div(vel))*v*dx
        elif method == 'gradq_avg':
            L2 = dot(avg(vel), n('+'))*jump(v)*dS - dot(vel, grad(v))*dx
        else:
            raise ValueError('Divergence type %r not supported' % method)
        
        # Pre-factorize matrices and store for usage in projection
        self._div_dS_solver = df.LocalSolver(a1, L1)
        self._div_dx_solver = df.LocalSolver(a2, L2)
        self._div_dS_solver.factorize()
        self._div_dx_solver.factorize()
        self._div_dS = df.Function(V)
        self._div_dx = df.Function(V)
        
    def _setup_energy(self, rho, vel, gvec, x0):
        """
        Calculate kinetic and potential energy
        """
        self._form_E_k = 1/2*rho*dot(vel, vel)*dx(domain=self.mesh)
        self._form_E_p = rho*dot(-gvec, x0)*dx(domain=self.mesh)
        
    def _setup_mass(self, rho):
        """
        Calculate mass
        """
        self._form_mass = rho*dx(domain=self.mesh)
    
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
    def divergences(self):
        """
        Calculate the difference between the flux on the same
        facet between two different cells and the divergence
        inside each cell.
        
        Returns the sum of facet errors for each cell and the
        divergence error in each cell as DG0 functions 
        """
        self._div_dS_solver.solve_global_rhs(self._div_dS)
        self._div_dx_solver.solve_global_rhs(self._div_dx)
        return self._div_dS, self._div_dx
    
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
