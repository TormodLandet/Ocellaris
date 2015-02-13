# encoding: utf8
"""
Solve the lid driven cavity flow using p1 DG elements
"""
from __future__ import division
import dolfin
from dolfin import dot, nabla_grad, avg, jump, dx, dS

###########################################################################################

dolfin.set_log_level(dolfin.WARNING)

# Mesh size and polynomial degree
N = 32
P = 1

# Time stepping
dt = 1/(8*N)
tmax = 1.0
num_inner_iter = 2
df_dt = dolfin.Constant(dt)

Re = 80
tramp = 0.2
nu1 = 1/Re # Start with this viscosity
nu2 = 1/Re # viscosity = nu2 after tramp seconds
umax = 1.0 # Velocity at lid 

# Physical properties used in simulation
rho = dolfin.Constant(1)
nu = dolfin.Constant(nu1)
g = dolfin.Constant([0, 0])
u_lid = dolfin.Constant(umax)

# Pure upwind differencing (beta=0 is upwind, beta=1 is downwind)
beta = dolfin.Constant(0.0)

# Geometry and function spaces
mesh = dolfin.UnitSquareMesh(N, N)
ndim = 2
Vu = dolfin.FunctionSpace(mesh, 'DG', P)
Vp = dolfin.FunctionSpace(mesh, 'DG', P)

# Show dolfin 3D plots?
PLOT = False

###########################################################################################

def define_advection_problem(V, up, upp, u_conv, n, beta, time_coeffs, dt, penalty, dirichlet_bcs):
    """
    Define the advection problem
    
     d/dt(u) + u_conv ⋅ grad(u) = 0
     
    Returns the bilinear and linear forms
    """
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
    c1, c2, c3 = time_coeffs[0],time_coeffs[1], time_coeffs[2] 
    eq = (c1*u + c2*up + c3*upp)/dt*v*dx \
         - u*dot(u_conv, nabla_grad(v))*dx \
         + flux*jump(v)*dS
    #     + u*dot(u_conv, n)*v*ds
    
    # Enforce Dirichlet BCs weakly
    for dbc in dirichlet_bcs:
        #eq += penalty*dot(u_conv, n)*v*(u - dbc.value)*dbc.ds
        eq += dot(u_conv, n)*v*dbc.value*dbc.ds
    
    a, L = dolfin.lhs(eq), dolfin.rhs(eq)
    
    return a, L

def define_poisson_problem(V, k, f, n, penalty, dirichlet_bcs, neumann_bcs):
    """
    Define the Poisson problem for u in f.space V
    
        - div(k*grad(u)) = f
    
    Note the minus in front of the first term!
    
    Returns the bilinear and linear forms
    """
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    
    # Interior
    a = dot(nabla_grad(v), nabla_grad(u))*dx
    
    # Interior facets
    a += avg(penalty)*jump(u)*jump(v)*dS \
         - dot(avg(nabla_grad(v)), n('+'))*jump(u)*dS \
         - dot(avg(nabla_grad(u)), n('+'))*jump(v)*dS

    # Source term
    L = -v*f/k*dx
    
    # Enforce Dirichlet BCs weakly
    for dbc in dirichlet_bcs:
        # These terms penalises jumps in (u - uD) on ds
        # where uD is the Dirichlet value on the boundary
        a += penalty*u*v*dbc.ds \
             - dot(nabla_grad(v), n)*u*dbc.ds \
             - dot(nabla_grad(u), n)*v*dbc.ds
        L += penalty*v*dbc.value*dbc.ds \
             - dot(nabla_grad(v), n)*dbc.value*dbc.ds
    
    # Enforce Neumann BCs weakly
    for nbc in neumann_bcs:
        L += v*nbc.value*nbc.ds
        
    # FIXME: introduce k in the equation properly!!!!
    a = k*a
    L = k*L
    
    return a, L

###########################################################################################

# Create velocity functions
u_list = [dolfin.Function(Vu) for _ in range(ndim)]
up_list = [dolfin.Function(Vu) for _ in range(ndim)]
upp_list = [dolfin.Function(Vu) for _ in range(ndim)]
u_conv_list = [dolfin.Function(Vu) for _ in range(ndim)]
u_star_list = [dolfin.Function(Vu) for _ in range(ndim)]
u_vec = dolfin.as_vector(u_list)
u_conv_vec = dolfin.as_vector(u_conv_list)
u_star_vec = dolfin.as_vector(u_star_list)

# Create pressure functions
p = dolfin.Function(Vp)
p_hat = dolfin.Function(Vp)

# Define boundary areas
class Walls(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
class Lid(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 0.999999
marker = dolfin.FacetFunction("size_t", mesh)
Walls().mark(marker, 1)
Lid().mark(marker, 2)

# Create a boundary measure that is aware of the marked regions
ds = dolfin.Measure('ds')[marker]

# Boundary conditions for the velocity
class DirichletBC(dolfin.DirichletBC):
    def __init__(self, V, value, subdomain_marker, subdomain_id, ds):
        super(DirichletBC, self).__init__(V, value, subdomain_marker, subdomain_id, method='geometric')
        self.ds = ds(subdomain_id)
        self.value = value
bc_walls = DirichletBC(Vu, dolfin.Constant(0), marker, 1, ds)
bc_lid_x = DirichletBC(Vu, u_lid, marker, 2, ds)
bc_lid_y = DirichletBC(Vu, dolfin.Constant(0), marker, 2, ds)
dirichlet_bcs_vel = [[bc_lid_x, bc_walls],
                     [bc_lid_y, bc_walls]]
neumann_bcs_vel = [[], []]

# Boundary conditions for the pressure 
class NeumannBC(object):
    def __init__(self, value, ds):
        self.value = value
        self.ds = ds
neumann_bcs_pres = [NeumannBC(dolfin.Constant(0), ds(1)),
                    NeumannBC(dolfin.Constant(0), ds(2))]
dirichlet_bcs_pres = []

# Mesh parameters
n = dolfin.FacetNormal(mesh)
h = dolfin.CellSize(mesh)

# Coefficients for u, up and upp 
time_coeffs = dolfin.Constant([1, -1, 0]) # First time step
time_coeffs2 = dolfin.Constant([3/2, -2, 1/2]) # All later time steps

# Define the momentum prediction equations
penalty = dolfin.Constant(100)/h
eqs_mom_pred = []
for d in range(ndim):
    f = -1/rho*p.dx(d) + g[d]
    a1, L1 = define_advection_problem(Vu, up_list[d], upp_list[d], u_conv_vec, n, beta, time_coeffs, df_dt, penalty, dirichlet_bcs_vel[d])
    a2, L2 = define_poisson_problem(Vu, nu, f, n, penalty, dirichlet_bcs_vel[d], neumann_bcs_vel[d])
    eq = a1+a2, L1+L2
    eqs_mom_pred.append(eq)

# Define the pressure correction equation
penalty = dolfin.Constant(100)/h
f = time_coeffs[0]/df_dt*dolfin.nabla_div(u_star_vec)
eq_pressure = define_poisson_problem(Vp, 1/rho, -f, n, penalty, dirichlet_bcs_pres, neumann_bcs_pres)

###########################################################################################

# Timer function
# A small decorator that stores the cummulative time spent in each 
# function that is wrapped by the decorator.
# Functions are identified by their names
from functools import wraps
import time
from collections import defaultdict
def timeit(f):
    """
    Decorator to time a set of functions
    and store the cummulative time
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        t1 = time.time()
        ret =  f(*args, **kwds)
        timeit.timings[f.__name__].append(time.time()-t1)
        return ret
    return wrapper
timeit.timings = defaultdict(list)

@timeit
def update_convection(t, dt):
    """
    Update explicit convective velocity used to linearise the convective term
    """
    for d in range(ndim):
        u_conv = u_conv_list[d]
        up =  up_list[d]
        upp = upp_list[d]
        u_conv.assign(2*up - upp)
        
        #for bc in dirichlet_bcs_vel[d]:
        #    bc.apply(uic.vector())

@timeit
def momentum_prediction(t, dt):
    """
    Solve the momentum prediction equation
    """
    for d in range(ndim):
        us = u_star_list[d]
        a, L = eqs_mom_pred[d]
        dolfin.solve(a == L, us) #, dirichlet_bcs_vel[d])

@timeit
def pressure_correction():
    """
    Solve the pressure correction equation
    
    We handle the case where only Neumann conditions are given
    for the pressure by taking out the nullspace, a constant shift
    of the pressure, by providing the nullspace to the solver
    """
    a, L = eq_pressure
    
    # Assemble matrix and vector
    A = dolfin.assemble(a)
    b = dolfin.assemble(L)
    
    # Create vector that spans the null space
    null_vec = dolfin.Vector(p_hat.vector())
    null_vec[:] = 1.0
    null_vec *= 1.0/null_vec.norm("l2")
    
    # Create null space basis object and attach to Krylov solver
    solver = dolfin.KrylovSolver(A, "gmres")
    null_space = dolfin.VectorSpaceBasis([null_vec])
    solver.set_nullspace(null_space)
    
    # Orthogonalize b with respect to the null space
    null_space.orthogonalize(b)
    
    solver.solve(p_hat.vector(), b)

@timeit
def velocity_update():
    """
    Update the velocity estimates with the updated pressure
    field from the pressure correction equation
    """
    c1 = time_coeffs[0]
    for d in range(ndim):
        us = u_star_list[d]
        u_new = u_list[d]
        
        # Update the velocity
        f = us - df_dt/(c1*rho) * p_hat.dx(d)
        un = dolfin.project(f, Vu)
        u_new.assign(un)
        
        #for bc in dirichlet_bcs_vel[d]:
        #    bc.apply(u_new.vector())

@timeit
def pressure_update():
    """
    Update the pressure
    """
    p.vector()[:] = p.vector()[:] + p_hat.vector()[:]
    
@timeit
def velocity_update_final():
    """
    Update the velocities at the end of the time step
    """
    for d in range(ndim):
        u = u_list[d]
        up = up_list[d]
        upp = upp_list[d]
        upp.assign(up)
        up.assign(u)

@timeit        
def errornorm(u1, u2):
    """
    The L-infinity error norm 
    """
    ua1 = u1.vector().array()
    ua2 = u2.vector().array()
    return abs(ua1 - ua2).max()

@timeit
def plotit():
    """
    Plot during simulations. Variables are stored on the function
    to make quick and dirty static variables to hold the plot
    information in between calls 
    """
    import numpy
    import matplotlib.pyplot as plt
    
    Nplots = 4
    if not hasattr(plotit, 'figs'):
        # This is the first call, lets set up the plots    
        plt.ion()
        plotit.figs = [plt.figure() for _ in range(Nplots)]
        plotit.axes = [fig.add_subplot(111) for fig in plotit.figs]
        plotit.lines = [ax.plot([], [], label='Re=%.0f'  % Re)[0] for ax in plotit.axes]
        
        # Data from Ghia, Ghia and Shin (1982)
        # Re = 100
        ghia_y = [1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5,
                  0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0]
        ghia_u = [1.0, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641,
                  -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.0]
        ghia_x = [1.0, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5,
                  0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0]
        ghia_v = [0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533,
                  0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.0]
        plotit.axes[0].plot(ghia_u, ghia_y, 'ks', label='Ghia Re=100')
        plotit.axes[0].legend(loc='lower right')
        plotit.axes[1].plot(ghia_x, ghia_v, 'ks', label='Ghia Re=100')
        plotit.axes[1].legend(loc='lower right')
    
    posvec = numpy.linspace(0, 1, 100)
    res = numpy.array([0.0])
    pos = numpy.array([0.0, 0.0])
    velx = numpy.zeros_like(posvec)
    vely = numpy.zeros_like(posvec)
    pres = numpy.zeros_like(posvec)
    for i, posi in enumerate(posvec):
        # Center plane in y-direction
        pos[:] = (0.5, posi)
        u_list[0].eval(res, pos)
        velx[i] = res[0]
        p.eval(res, pos)
        pres[i] = res[0]
        # Center plane in x-direction
        pos[:] = (posi, 0.5)
        u_list[1].eval(res, pos)
        vely[i] = res[0]
    
    Q = 0
    plotit.lines[Q].set_xdata(velx)
    plotit.lines[Q].set_ydata(posvec)
    plotit.axes[Q].set_xlim(-0.5, 1)
    plotit.axes[Q].set_ylim(0, 1)
    plotit.axes[Q].set_title('Velocity in x at x=0.5')
    plotit.axes[Q].set_xlabel('U [m/s]')
    plotit.axes[Q].set_ylabel('y [m]')
    
    Q += 1
    plotit.lines[Q].set_ydata(vely)
    plotit.lines[Q].set_xdata(posvec)
    plotit.axes[Q].set_ylim(-0.5, 0.5)
    plotit.axes[Q].set_xlim(0, 1)
    plotit.axes[Q].set_title('Velocity in y at y=0.5')
    plotit.axes[Q].set_ylabel('U [m/s]')
    plotit.axes[Q].set_xlabel('y [m]')
    
    Q += 1
    plotit.lines[Q].set_xdata(pres)
    plotit.lines[Q].set_ydata(posvec)
    plotit.axes[Q].set_xlim(-0.6, 0.6)
    plotit.axes[Q].set_ylim(0, 1)
    plotit.axes[Q].set_title('Pressure at x=0.5')
    plotit.axes[Q].set_xlabel('P [Pa]')
    plotit.axes[Q].set_ylabel('y [m]')
    
    Q += 1
    tlist = list(plotit.lines[Q].get_xdata()) + [t]
    ulist = list(plotit.lines[Q].get_ydata()) + [maxvel]
    plotit.lines[Q].set_xdata(tlist)
    plotit.lines[Q].set_ydata(ulist)
    plotit.axes[Q].relim()
    plotit.axes[Q].autoscale_view()
    plotit.axes[Q].set_xlim(0, tmax)
    plotit.axes[Q].set_title('Maximum velocity magnitude')
    plotit.axes[Q].set_xlabel('t [s]')
    plotit.axes[Q].set_ylabel('U [m/s]')
    
    for fig in plotit.figs:
        fig.canvas.draw()
        fig.canvas.flush_events()

###########################################################################################

# Run the simulation
t = 0
it = 0
while t+dt <= tmax + 1e-8:
    it += 1
    t += dt
    print 'TIMESTEP %5d  -  t = %10.4f' % (it, t)
    
    #g.assign(dolfin.Constant([0, -min(1, t/tmax*2)]))
    #u_lid.assign(umax*min(1, t/tramp))
    if t > tramp:
        nu.assign(nu2)
    
    update_convection(t, dt)
    
    for inner_iter in range(num_inner_iter):
        momentum_prediction(t, dt)
        pressure_correction()
        velocity_update()
        pressure_update()
        
        # DEBUG:
        u0a = u_list[0].vector().array()
        u1a = u_list[1].vector().array()
        maxvel = (u0a**2+u1a**2).max()**0.5
        print '        Inner iter %d, max velocity %15.4e' % (inner_iter, maxvel),
        error = sum(errornorm(u_list[d], u_star_list[d]) for d in range(ndim))
        print ' errornorm*: %10.3e' % error,
        errorC = sum(errornorm(u_list[d], u_conv_list[d]) for d in range(ndim))
        print ' errornormC: %10.3e' % errorC
    
    velocity_update_final()
    plotit()
    
    # Change to the second order time derivative for velocities
    time_coeffs.assign(time_coeffs2)
    
    # DEBUG:
    pa = p.vector().array()
    print '    Velocity in x direction ranging from', u0a.min(), 'to', u0a.max()
    print '    Velocity in y direction ranging from', u1a.min(), 'to', u1a.max()
    print '    Pressure ranging from', pa.min(), 'to', pa.max()
    if maxvel > 1.5:
        print 'Diverging ...'
        break

# Print the runtime of the functions timed with the @timeit decorator
print '\n\nTIMINGS:'
tottime = sum(sum(vals) for vals in timeit.timings.itervalues())
for funcname in sorted(timeit.timings):
    durations = timeit.timings[funcname]
    print '%30s total time %7.3fs for %5d runs, minimum runtime %7.3fs' % \
          (funcname, sum(durations), len(durations), min(durations)),
    print '  (%5.1f%% of tot.time)' % (sum(durations)/tottime*100) 

###########################################################################################

if PLOT:
    dolfin.plot(marker, title='boundary')
    dolfin.plot(u_vec, title='u')
    dolfin.plot(u_list[0], title='u_x')
    dolfin.plot(u_list[1], title='u_y')
    dolfin.plot(p, title='p')
    dolfin.interactive()

raw_input('Press any key to quit')

