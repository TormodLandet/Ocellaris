from dolfin import Function, LocalSolver


class RungeKuttaDGTimestepping(object):
    def __init__(self, simulation, a, L, u, up, order=None):
        """
        RKDG timestepping. A is a block diagonal mass matrix form (u*v*dx),
        L is the form of the right hand side of du/dt = L. The functions
        u and up are the result and the previous value respectively. If order
        is not given it is taken as P+1 where P is the element degree of u.
        
        It is assumed that up is used in the form L so that changing it causes
        L to change. L should not depend on u or t, only up. The up function
        *will be modified* when running .step()! It will not contain a sensible
        value afterwords, while the u function will contain the time integrated
        value after the time step. The up function should contain the value of
        u at the beginning of the time step before running .step()
        """
        self.simulation = simulation
        
        V = u.function_space()
        if order is None:
            order = V.ufl_element().degree() + 0
        self.order = order
        
        # Number of stages
        if order <= 3:
            S = order
        else:
            self.A, self.B, _C = get_ssp_rk_coefficients(order)
            S = len(self.A)
        simulation.log.info('    Preparing SSP RK method of order %d with %d stages' % (order, S))
        
        self.u = u
        self.up = up
        self.us = [Function(V) for _ in range(S-1)]
        self.dus = [Function(V) for _ in range(S)]
        self.solver = LocalSolver(a, L)
        self.solver.factorize()
        
        if self.order > 3:
            # Need one extra function storage when running the generic code
            self.up_store = Function(V)
            
    def _solve(self, du, fdt, uexpl):
        if uexpl is not self.up:
            self.up.assign(uexpl)
        self.solver.solve_local_rhs(du)
    
    def step(self, dt):
        """
        Use Runge-Kutta to step dt forward in time
        """
        u, up = self.u, self.up
        us, dus = self.us, self.dus
        K, ls = self.order, self.solver
        
        if K == 1:
            u.assign(up)
            ls.solve_global_rhs(dus[0])
            u.vector().axpy(dt, dus[0].vector())
        
        elif K == 2:
            u.assign(up)
            self._solve(dus[0], 0.0, up)
            us[0].assign(u)
            us[0].vector().axpy(dt, dus[0].vector())
            
            self._solve(dus[1], 1.0, us[0])
            u.vector().axpy(0.5*dt, dus[0].vector())
            u.vector().axpy(0.5*dt, dus[1].vector())
        
        elif K == 3:
            u.assign(up)
            self._solve(dus[0], 0.0, up)
            us[0].assign(u)
            us[0].vector().axpy(dt, dus[0].vector())
            
            self._solve(dus[1], 0.5, us[0])
            us[1].assign(u)
            us[1].vector().axpy(0.25*dt, dus[0].vector())
            us[1].vector().axpy(0.25*dt, dus[1].vector())
            
            self._solve(dus[2], 1.0, us[1])
            u.vector().axpy(1/6*dt, dus[0].vector())
            u.vector().axpy(1/6*dt, dus[1].vector())
            u.vector().axpy(2/3*dt, dus[2].vector())
            
        else:
            # Generic implementation
            self.up_store.assign(up)
            A, B = self.A, self.B
            S = len(A) # number of stages
            
            for i in range(S):
                if i > 0:
                    up.assign(us[i-1])
                ls.solve_local_rhs(dus[i])

                un = us[i] if i < S-1 else u
                un.vector().zero()

                for j in range(i+1):        
                    a, b = A[i][j], B[i][j]

                    if a != 0:
                        uo = us[j-1] if j > 0 else self.up_store
                        un.vector().axpy(a, uo.vector())
                    if b != 0:
                        un.vector().axpy(b*dt, dus[j].vector())


def get_ssp_rk_coefficients(K):
    """
    Get coefficients for the strong stability-preserving
    Runge-Kutta method of order K. May return more stages
    S than the order K to avid negative B coeffients and
    provide more CFL stability (S >= K).
    
    Returns RK SSP coefficients A, B and the effective CFL
    multiplier C
    
    See i.e
    Ruuth - 2005 - "Global optimization of Explicit
      Strong-Stability-Preserving Runge-Kutta methods",
    Kubatko et al - 2008 - "Time step restriction for
      Runge-Kutta discontinuous Galerkin methods on
      triangular grids",
    and several others starting with Shu & Osher (1988)
    """
    C = 1.0
    if K == 1:
        # SSP(1,1): Forward Euler
        A = [[1]]
        B = [[1]]
    elif K == 2:
        # SSP(2,2): Two stage RK2 scheme
        A = [[1], [1/2, 1/2]]
        B = [[1], [0, 1/2]]
    elif K == 3:
        # SSP(3,3): Three stage RK3 SSP scheme
        A = [[1], [3/4, 1/4], [1/3, 0, 2/3]]
        B = [[1], [0, 1/4], [0, 0, 2/3]]
    elif K == 4:
        # SSP(5,4): five stage RK4 SSP scheme
        A = [[1.0], [0.444370493651235, 0.555629506348765],
            [0.620101851488403, 0, 0.379898148511597],
            [0.178079954393132, 0, 0, 0.821920045606868],
            [0, 0, 0.517231671970585, 0.096059710526147, 0.386708617503269]]
        B = [[0.391752226571890], [0, 0.368410593050371],
            [0, 0, 0.251891774271694], [0, 0, 0, 0.544974750228521],
            [0, 0, 0, 0.063692468666290, 0.226007483236906]]
        C = 1.50818004918983
    else:
        raise NotImplementedError('No SSP RK method of order %d implemented' % K)
    return A, B, C
