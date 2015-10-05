Convergence rates
=================

Taylor-Green
------------

The `Taylor-Green vortex <http://en.wikipedia.org/wiki/Taylor-Green_vortex>`_
is an analytical solution of the incompressible Navier-Stokes equations in a
doubly periodic 2D domain with side lengths equal to 2.0. The solution is given
for :math:`t \ge 0` as:

.. math::

    u_0 &= -\sin(\pi x_1)\cos(\pi x_0) e^{-2\pi^2 \nu t} \\
    u_1 &= \quad\sin(\pi x_0)\cos(\pi x_1) e^{-2\pi^2 \nu t} \\
    p   &= -\frac{1}{4}(\cos(2\pi x_0) + \cos(2 \pi x_1))  e^{-4 \pi^2\nu t}

To study the convergence order of the Ocellaris solver we calculate the
:math:`L_2` norm of the difference between the analytical solution and the
calculated solutions after 1 second simulation and with kinematic viscosity
:math:`\nu=0.01`.


Convergence in space
....................

We study the order of convergence in space which we expect to be 
:math:`\mathcal{O}(N^{-k})` where :math:`N` is the number of elements along
each side of the mesh and :math:`k` is the convergence rate. We keep the time
step constant :math:`\Delta t=0.001` and vary the number of elements in each 
direction. The convergence rate is calulated by comparing two consecutive
discretisations and taking

.. math::

    k = \frac {||u-u_h||^{N_i}}{||u-u_h||^{N_{i+1}}} \log(\frac{N_{i+1}}{N_i})^{-1} 

To test the continuous Galerkin FEM implementation of the IPCS solver we use
P2P1 Tayloor-Hood elements for the spatial discretisation. This means that the
velocity is approximated with quadratic polynomials while the pressure is 
approximated with linear polynomials.

.. table:: Results with CG P2P1 elements, IPCS solver with BDF timestepping

    ===== =============== ==============  ===== =====
        N      L2 error u     L2 error p    k u   k p
    ===== =============== ==============  ===== =====
        8       8.881e-02      1.544e-01      -     -
       16       5.202e-03      4.973e-02   4.09  1.63
       24       9.712e-04      2.278e-02   4.14  1.93
       32       3.008e-04      1.286e-02   4.07  1.99
       40       1.219e-04      8.233e-03   4.05  2.00
    ===== =============== ==============  ===== =====
    
The convergence rate in the velocities :math:`k=4` can not be expected on a
general mesh. In fact, if we take the same mesh vertices and use alternating
diagonals instead of only right pointing diagonals in the mesh we get a
convergence rate which is just a bit higher than the expected rate :math:`k=3`:

.. table:: Results with CG P2P1 elements, IPCS solver with BDF timestepping
    
    ===== =============== ===============  ===== =====
        N      L2 error u      L2 error p   k u   k p
    ===== =============== ===============  ===== =====
        8       8.802e-02       1.696e-01      -     -
       16       6.637e-03       5.008e-02   3.73  1.76
       24       1.558e-03       2.290e-02  3.57  1.93
       32       5.965e-04       1.295e-02   3.34  1.98
       40       2.915e-04       8.297e-03   3.21  1.99
    ===== =============== ===============  ===== =====

When using discontinuous Galerkin elements we need to work around the fact that
FEniCS 1.5 does not support periodic boundary conditions for the DG function 
space. We do this by imposing Dirichlet conditions on all fields on all
boundaries. The results are as follows when not using the superconvergent mesh:

.. table:: Results with DG P2P1 elements, IPCS solver with BDF timestepping. Dirichlet BCs
    
    ===== =============== ===============  ===== =====
        N      L2 error u      L2 error p    k u   k p
    ===== =============== ===============  ===== =====
        8       7.808e-02       2.755e-01      -     -
       16       7.814e-03       6.099e-02   3.32  2.18
       24       1.830e-03       2.953e-02   3.58  1.79
       32       6.970e-04       1.719e-02   3.35  1.88
       40       3.411e-04       1.116e-02   3.20  1.94
    ===== =============== ===============  ===== =====
    
To verify that the use of Dirichlet boundary conditions does not significantly
skew the results we also show the results from using the continuous elements
with pure Dirichlet boundary conditions:

.. table:: Results with CG P2P1 elements, IPCS solver with BDF timestepping. Dirichlet BCs
 
    =====  =============== ===============  ===== =====
        N      L2 error u1      L2 error p   k u1   k p
    =====  =============== ===============  ===== =====
        8        1.050e-01       1.504e-01      -     -
       16        1.423e-02       4.860e-02   2.88  1.63
       24        4.788e-03       2.446e-02   2.69  1.69
       32        2.270e-03       1.452e-02   2.60  1.81
       40        1.269e-03       9.620e-03   2.61  1.85
    =====  =============== ===============  ===== =====

.. warning:: Find out why this happened


Convergence in time
...................

To study convergence in time it is crucial to eliminate the errors in the
spatial discretisation. Here this is done by having N=60 and using P4P3 
elements (polynomial degree 4 for the velocities and 3 for the pressure). This
ensures that we are only studing the error in the time discretisation.

The temporal convergence is assumed to be :math:`\mathcal{O}(\Delta t^k)` where
:math:`\Delta t` is the time step and :math:`k` is the convergence rate.

.. table:: Results with CG P4P3 elements, IPCS solver with BDF timestepping

    ======= =============== ==============  ===== =====
         dt      L2 error u     L2 error p    k u   k p
    ======= =============== ==============  ===== =====
    0.50000       4.447e-04      1.166e-02      -     -
    0.25000       1.441e-04      2.850e-03   1.63  2.03
    0.12500       3.815e-05      7.009e-04   1.92  2.02
    0.06250       9.786e-06      1.740e-04   1.96  2.01
    0.03120       2.467e-06      4.339e-05   1.98  2.00
    ======= =============== ==============  ===== =====


.. _sec-taylor-green-convergence-script:

Script to reproduce results
...........................

The convergence rates above where calculated with the following script:

.. code-block:: python
    
    from __future__ import division
    import time
    from math import log
    import dolfin
    from ocellaris import Simulation, run_simulation
    
    def run_and_calculate_error(N, dt, tmax, polydeg_u, polydeg_p, num_inner_iter, timestepping_method):
        """
        Run Ocellaris and return L2 errors in the last time step
        """
        # Setup and run simulation
        sim = Simulation()
        sim.input.read_yaml('taylor-green.inp')
        sim.input['mesh']['Nx'] = sim.input['mesh']['Ny'] = N
        sim.input['mesh']['diagonal'] = 'right' # Gives superconvergence in space
        #sim.input['mesh']['diagonal'] = 'right/left' # Gives closer to normal convergence in space
        sim.input['time']['dt'] = dt
        sim.input['time']['tmax'] = tmax
        sim.input['solver']['polynomial_degree_velocity'] = polydeg_u
        sim.input['solver']['polynomial_degree_pressure'] = polydeg_p
        sim.input['solver']['timestepping_method'] = timestepping_method
        sim.input['solver']['num_inner_iter'] = num_inner_iter
        run_simulation(sim)
        
        # Interpolate the analytical solution to the same function space
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']
        vals = dict(t=sim.time, nu=sim.input['physical_properties']['nu0'])
        u0e = dolfin.Expression('-sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)', **vals)
        u1e = dolfin.Expression(' sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)', **vals)
        pe  = dolfin.Expression('-(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*t)/4', **vals)
        u0a = dolfin.interpolate(u0e, Vu)
        u1a = dolfin.interpolate(u1e, Vu)
        pa = dolfin.interpolate(pe, Vp)
        
        # Calculate errors
        err_u0 = calc_err(sim.data['u0'], u0a)
        err_u1 = calc_err(sim.data['u1'], u1a)
        err_p = calc_err(sim.data['p'], pa)
        
        return err_u0, err_u1, err_p
    
    def calc_err(f_num, f_ana):
        """
        Calculate scaled L2 error
        """
        f_err = dolfin.Function(f_num.function_space())
        f_err.vector()[:] = f_ana.vector()[:] - f_num.vector()[:]
        return dolfin.norm(f_err) / dolfin.norm(f_ana) 
    
    def run_convergence_space(N_list):
        dt = 0.001
        tmax = 1.0
        num_inner_iter = 3
        results = {}
        for N in N_list:
            t1 = time.time()
            print 'Running N =', N
            err_u0, err_u1, err_p = run_and_calculate_error(N=N, dt=dt, tmax=tmax,
                                                            polydeg_u=2, polydeg_p=1,
                                                            num_inner_iter=num_inner_iter,
                                                            timestepping_method='BDF')
            results[N] = err_u0, err_u1, err_p
            
            print 'Time spent:', time.time() - t1
            print 'err_u0_max %15.4e' % err_u0
            print 'err_u1_max %15.4e' % err_u1
            print 'err_p_max  %15.4e' % err_p
            print
            
        print 'dt =', dt
        print '===== =============== =============== ===============  ===== ===== ====='
        print '    N     L2 error u0     L2 error u1      L2 error p   k u0  k u1   k p'
        print '===== =============== =============== ===============  ===== ===== ====='
        for i, N in enumerate(N_list):
            eu0, eu1, ep = results[N]
            if i == 0:
                print '%5d %15.3e %15.3e %15.3e' % (N, eu0, eu1, ep)
            else:
                prev_N = N_list[i-1]
                prev_eu0, prev_eu1, prev_ep = results[prev_N]
                fac = log(N/prev_N)
                print '%5d %15.3e %15.3e %15.3e ' % (N, eu0, eu1, ep),
                print '%5.2f %5.2f %5.2f' % (log(prev_eu0/eu0)/fac, log(prev_eu1/eu1)/fac, log(prev_ep/ep)/fac)
        print '===== =============== =============== ===============  ===== ===== ====='        
    
    
    def run_convergence_time(dt_list):
        N = 60
        tmax = 1.0
        num_inner_iter = 20
        results = {}
        for dt in dt_list:
            t1 = time.time()
            print 'Running dt =', dt
            err_u0, err_u1, err_p = run_and_calculate_error(N=N, dt=dt, tmax=tmax,
                                                            polydeg_u=4, polydeg_p=3,
                                                            num_inner_iter=num_inner_iter,
                                                            timestepping_method='BDF')
            results[dt] = err_u0, err_u1, err_p
            
            print 'Time spent:', time.time() - t1
            print 'err_u0_max %15.4e' % err_u0
            print 'err_u1_max %15.4e' % err_u1
            print 'err_p_max  %15.4e' % err_p
            print
            
        print 'N =', N
        print '======= =============== =============== ===============  ===== ===== ====='
        print '     dt     L2 error u0     L2 error u1      L2 error p   k u0  k u1   k p'
        print '======= =============== =============== ===============  ===== ===== ====='
        for i, dt in enumerate(dt_list):
            eu0, eu1, ep = results[dt]
            if i == 0:
                print '%7.5f %15.3e %15.3e %15.3e' % (dt, eu0, eu1, ep)
            else:
                prev_dt = dt_list[i-1]
                prev_eu0, prev_eu1, prev_ep = results[prev_dt]
                fac = log(prev_dt/dt)
                print '%7.5f %15.3e %15.3e %15.3e ' % (dt, eu0, eu1, ep),
                print '%5.2f %5.2f %5.2f' % (log(prev_eu0/eu0)/fac, log(prev_eu1/eu1)/fac, log(prev_ep/ep)/fac)
        print '======= =============== =============== ===============  ===== ===== ====='        
    
    
    run_convergence_time([5e-1, 2.5e-1, 1.25e-1, 6.25e-2, 3.12e-2])
    run_convergence_space([8, 16, 24, 32, 40])
    
The input file ``taylor-green.inp`` looks as follows. Note that the boundary
conditions must be changed from periodic to Dirichlet in order to use DG
elements. This is done by commenting the periodic boundary conditions and 
uncommenting the Dirichlet boundary conditions in the file below. 

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0
        
    metadata:
        author: Tormod Landet
        date: 2015-03-13
        description: |
            Implements the Taylor-Green vortex test case. This benchmark case
            with purely periodic boundary conditions has an analytical solution
            in both space and time with the incompressible Navier-Stokes equations
    
    physical_properties:
        g: [0, 0]
        nu0: 0.01
        rho0: 1.0
    
    mesh:
        type: Rectangle
        Nx: 64
        Ny: 64
        endx: 2
        endy: 2
    
    # Periodic boundary conditions. Not supported in FEniCS 1.5 for DG. Works for CG
    boundary_conditions:
    -   name: left and bottom    
        selector: code
        inside_code: |
            inside = bool((near(x[0], 0) or near(x[1], 0)) and 
                          (not ((near(x[0], 0) and near(x[1], 2)) or 
                          (near(x[0], 2) and near(x[1], 0)))) and on_boundary)
        map_code: |
            if near(x[0], 2) and near(x[1], 2):
                y[0] = x[0] - 2.0
                y[1] = x[1] - 2.0
            elif near(x[0], 2):
                y[0] = x[0] - 2.0
                y[1] = x[1]
            else:
                y[0] = x[0]
                y[1] = x[1] - 2.0
    
    # Dirichlet boundary conditions for all variables
    #boundary_conditions:
    #-   name: walls
    #    selector: code
    #    inside_code: on_boundary
    #    u:
    #        type: CppCodedValue
    #        cpp_code:
    #        -   -sin(pi*x[1]) * cos(pi*x[0]) * exp(-2*pi*pi*nu*t)
    #        -    sin(pi*x[0]) * cos(pi*x[1]) * exp(-2*pi*pi*nu*t)
    #    p:
    #        type: CppCodedValue
    #        cpp_code: -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4.*pi*pi*nu*t)/4
    
    initial_conditions:
        up0:
            cpp_code: -sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)
        up1:
            cpp_code:  sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)
        upp0:
            cpp_code: -sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*(t-dt))
        upp1:
            cpp_code:  sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*(t-dt))
        p:
            cpp_code: -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*t)/4
    
    time:
        dt: 0.001
        tmax: 1.0
    
    output:
        prefix: taylor_green
        log_name: .log
        dolfin_log_level: warning
        ocellaris_log_level: warning
    
    solver:
        type: IPCS
        polynomial_degree_pressure: 1
        polynomial_degree_velocity: 2
        function_space_pressure: CG
        function_space_velocity: CG
        timestepping_method: BDF
        u:
            parameters:
                relative_tolerance: 1.0e-10
                absolute_tolerance: 1.0e-15
        p:
            parameters:
                relative_tolerance: 1.0e-10
                absolute_tolerance: 1.0e-15
 