# encoding: utf8
"""
Solve the lid driven cavity flow using p1 DG elements
"""
from __future__ import division
from functools import wraps
import time
from collections import defaultdict
import numpy
from matplotlib import pyplot
import dolfin
from dolfin import dot, nabla_grad, nabla_div, avg, jump, dx, dS
from ocellaris.utils import debug_console_hook

###########################################################################################

dolfin.set_log_level(dolfin.WARNING)

# Mesh size and polynomial degree
N = 64
Pu = 2
Pp = 1
rotational = False  

# Time stepping
dt = 0.05 # 1/(8*N)
tmax = 40.0
min_inner_iter = 1
max_inner_iter = 3
u_star_error_limit = 1e-2
df_dt = dolfin.Constant(dt)

Re = 1000
umax = 1.0 # Velocity at lid
kin_visc = 1/Re
max_allowed_velocity = 10

# Physical properties used in simulation
rho = dolfin.Constant(1)
nu = dolfin.Constant(kin_visc)
g = dolfin.Constant([0, 0])
u_lid = dolfin.Constant(umax)

# Pure upwind differencing (beta=0 is upwind, beta=1 is downwind)
beta = dolfin.Constant(0.0)

# Geometry and function spaces
mesh = dolfin.UnitSquareMesh(N, N)

# Grade mesh towards the walls
if False:
    x = mesh.coordinates()
    x[:] = (x - 0.5) * 2
    x[:] = 0.5*(numpy.cos(numpy.pi*(x-1.) / 2.) + 1.)

ndim = 2
Vu = dolfin.FunctionSpace(mesh, 'CG', Pu)
Vp = dolfin.FunctionSpace(mesh, 'CG', Pp)

# Show dolfin 3D plots?
PLOT = True

###########################################################################################

def define_advection_problem(u, v, up, upp, u_conv, n, beta, time_coeffs, dt, penalty, dirichlet_bcs):
    """
    Define the advection problem
    
     d/dt(u) + u_conv ⋅ grad(u) = 0
     
    Returns the bilinear and linear forms
    """
    family = u.element().family()
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        c1, c2, c3 = time_coeffs 
        eq = (c1*u + c2*up + c3*upp)/dt*v*dx + dot(u_conv, nabla_grad(u))*v*dx
        a, L = dolfin.lhs(eq), dolfin.rhs(eq)

    elif family == 'Discontinuous Lagrange':
        # Upstream and downstream normal velocities
        flux_nU = u*(dot(u_conv, n) + abs(dot(u_conv, n)))/2
        flux_nD = u*(dot(u_conv, n) - abs(dot(u_conv, n)))/2
        
        # Define the blended flux
        # The blending factor beta is not DG, so beta('+') == beta('-')
        b = beta('+')
        flux = (1-b)*(flux_nU('+') - flux_nU('-')) + b*(flux_nD('+') - flux_nD('-'))
        
        # Equation to solve
        c1, c2, c3 = time_coeffs 
        eq = (c1*u + c2*up + c3*upp)/dt*v*dx \
             - u*dot(u_conv, nabla_grad(v))*dx \
             + flux*jump(v)*dS
        
        # Enforce Dirichlet BCs weakly
        for dbc in dirichlet_bcs:
            #eq += penalty*dot(u_conv, n)*v*(u - dbc.value)*dbc.ds
            eq += dot(u_conv, n)*v*(u - dbc.value)*dbc.ds
        
        a, L = dolfin.lhs(eq), dolfin.rhs(eq)
        
    return a, L

def define_poisson_problem(u, v, k, f, n, penalty, dirichlet_bcs, neumann_bcs):
    """
    Define the Poisson problem for u in f.space V
    
        - div(k*grad(u)) = f
    
    Note the minus in front of the first term!
    
    Returns the bilinear and linear forms
    """
    family = u.element().family()
    
    if family == 'Lagrange':
        # Continous Galerkin implementation 
        a = k*dot(nabla_grad(v), nabla_grad(u))*dx
        L = v*f*dx
        # Enforce Neumann BCs weakly
        for nbc in neumann_bcs:
            L += k*v*nbc.value*nbc.ds
    
    elif family == 'Discontinuous Lagrange':
        # Discontinous Galerkin implementation (Symmetric Interior Penalty method)
        # See  Epshteyn and Rivière, 2007:
        # "Estimation of penalty parameters for symmetric interior penalty Galerkin methods"
        a = k*dot(nabla_grad(u), nabla_grad(v))*dx
        a -= avg(k*dot(n, nabla_grad(u)))*jump(v)*dS
        a -= avg(k*dot(n, nabla_grad(v)))*jump(u)*dS
        a += avg(penalty)*jump(u)*jump(v)*dS
        L = f*v*dx
        
        for dbc in dirichlet_bcs:
            a -= k*dot(n, nabla_grad(u))*v*dbc.ds
            a -= k*dot(n, nabla_grad(v))*u*dbc.ds
            a += penalty*u*v*dbc.ds
            L += dbc.value*(penalty*v - k*dot(nabla_grad(v), n))*dbc.ds
        
        for nbc in neumann_bcs:
            L += nbc.value*v*nbc.ds
    
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
        return on_boundary and x[1] > 0.999999 #and x[0]*(1-x[0]) > 0.00001
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
neumann_bcs_pres = [NeumannBC(dolfin.Constant(0), ds(1)), NeumannBC(dolfin.Constant(0), ds(2))]
dirichlet_bcs_pres = [] #DirichletBC(Vu, dolfin.Constant(0), marker, 2, ds)]

# Mesh parameters
n = dolfin.FacetNormal(mesh)
h = dolfin.CellSize(mesh)
vol = dolfin.CellVolume(mesh)
area = dolfin.FacetArea(mesh)

# Coefficients for u, up and upp 
time_coeffs = dolfin.Constant([1, -1, 0]) # First time step
time_coeffs2 = dolfin.Constant([3/2, -2, 1/2]) # All later time steps

# Define the momentum prediction equations
#penalty = dolfin.Constant(100)/h
#penalty = (Pu +1)*(Pu + ndim)/ndim*2*area*(ndim + 1)/vol

geom_fac = 0
for cell in dolfin.cells(mesh):
    vol = cell.volume()
    area = sum(cell.facet_area(i) for i in range(ndim+1))
    gf = area/vol
    geom_fac = max(geom_fac, gf)
geom_fac *= 1.0

penalty = dolfin.Constant((Pu + 1)*(Pu + ndim)/ndim * geom_fac)
eqs_mom_pred = []
for d in range(ndim):
    trial = dolfin.TrialFunction(Vu)
    test = dolfin.TestFunction(Vu)
    f = -1/rho*p.dx(d) + g[d]
    a1, L1 = define_advection_problem(trial, test, up_list[d], upp_list[d],
                                      u_conv_vec, n, beta, time_coeffs, df_dt,
                                      penalty, dirichlet_bcs_vel[d])
    a2, L2 = define_poisson_problem(trial, test, nu, f, n, penalty, dirichlet_bcs_vel[d], neumann_bcs_vel[d])
    eq = a1+a2, L1+L2
    eqs_mom_pred.append(eq)

# Define the pressure correction equation
#penalty = dolfin.Constant(10000)/h
penalty = dolfin.Constant((Pp + 1)*(Pp + ndim)/ndim * geom_fac)
trial = dolfin.TrialFunction(Vp)
test = dolfin.TestFunction(Vp)
f = time_coeffs[0]/df_dt * nabla_div(u_star_vec)
eq_pressure = define_poisson_problem(trial, test, 1/rho, -f, n, penalty, dirichlet_bcs_pres, neumann_bcs_pres)

###########################################################################################

# Timer function
# A small decorator that stores the cummulative time spent in each
# function that is wrapped by the decorator.
# Functions are identified by their names
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
        
        for bc in dirichlet_bcs_vel[d]:
            bc.apply(u_conv.vector())

@timeit
def momentum_prediction(t, dt):
    """
    Solve the momentum prediction equation
    """
    for d in range(ndim):
        us = u_star_list[d]
        a, L = eqs_mom_pred[d]
        if us.element().family() == 'Lagrange':
            dolfin.solve(a == L, us, dirichlet_bcs_vel[d])
        else:
            dolfin.solve(a == L, us)

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
    
    if p_hat.element().family() == 'Lagrange':
        for bc in dirichlet_bcs_pres:
            bc.apply(A, b)
    
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
        
        for bc in dirichlet_bcs_vel[d]:
            bc.apply(u_new.vector())

@timeit
def pressure_update():
    """
    Update the pressure
    """
    if rotational:
        f = p_hat - nu*nabla_div(u_star_vec)
        ph = dolfin.project(f, Vp)
    else:
        ph = p_hat
    p.vector()[:] = p.vector()[:] + ph.vector()[:]

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
    The L-infinity error norm of a vector valued function
    """
    res = 0
    for d in range(ndim):
        a1 = u1[d].vector().array()
        a2 = u2[d].vector().array()
        err = numpy.abs(a1 - a2).max()  
        res = max(res, err)
    return res

# Data from Ghia, Ghia and Shin (1982)
ghia_x = [1.0, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5,
                  0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0]
ghia_y = [1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5,
          0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0]
ghia_data = {100: 
                {'u': [1.0, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641,
                       -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.0],
                 'v': [0, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533,
                       0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.0]},
             1000:
                {'u': [1, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648,
                       -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0],
                 'v': [0, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 
                       0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.0]}}

@timeit
def ghia_error():
    """
    Calculate L2 error wrt Ghia data
    """
    if Re not in ghia_data:
        return -1
    
    
    error = 0
    res = numpy.array([0.0])
    pos = numpy.array([0.5, 0.0])
    for i, y in enumerate(ghia_y):
        pos[1] = y
        # Get calculated value
        u_list[0].eval(res, pos)
        calc_val = res[0]
        # Get Ghia et als value
        ghia_val = ghia_data[Re]['u'][i]
        # Sum squared errors 
        error += (ghia_val - calc_val)**2
    return error**0.5

@timeit
def plotit():
    """
    Plot during simulations. Variables are stored on the function
    to make quick and dirty static variables to hold the plot
    information in between calls 
    """
    Nplots = 2
    if not hasattr(plotit, 'figs'):
        # This is the first call, lets set up the plots    
        pyplot.ion()
        plotit.figs = [pyplot.figure() for _ in range(Nplots)]
        plotit.axes = [fig.add_subplot(111) for fig in plotit.figs]
        plotit.lines = [ax.plot([], [], label='Re=%.0f'  % Re)[0] for ax in plotit.axes]
        
        if Re in ghia_data:
            ghia_u = ghia_data[Re]['u']
            ghia_v = ghia_data[Re]['v']
            plotit.axes[0].plot(ghia_u, ghia_y, 'ks', label='Ghia Re=%d' % Re)
            plotit.axes[0].legend(loc='lower right')
            plotit.axes[1].plot(ghia_x, ghia_v, 'ks', label='Ghia Re=%d' % Re)
            plotit.axes[1].legend(loc='lower left')
    
    posvec = numpy.linspace(0, 1, 1000)
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
    if Q < Nplots:
        plotit.lines[Q].set_ydata(vely)
        plotit.lines[Q].set_xdata(posvec)
        plotit.axes[Q].set_ylim(-0.6, 0.6)
        plotit.axes[Q].set_xlim(0, 1)
        plotit.axes[Q].set_title('Velocity in y at y=0.5')
        plotit.axes[Q].set_ylabel('U [m/s]')
        plotit.axes[Q].set_xlabel('y [m]')
    
    Q += 1
    if Q < Nplots:
        plotit.lines[Q].set_xdata(pres)
        plotit.lines[Q].set_ydata(posvec)
        plotit.axes[Q].set_xlim(-0.6, 0.6)
        plotit.axes[Q].set_ylim(0, 1)
        plotit.axes[Q].set_title('Pressure at x=0.5')
        plotit.axes[Q].set_xlabel('P [Pa]')
        plotit.axes[Q].set_ylabel('y [m]')
    
    Q += 1
    if Q < Nplots:
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

def debug_console_hook_wrapper():
    class Dummy():
        pass
    simulation = Dummy()
    simulation.data = {k:v for k, v in globals().items() if k not in dolfin.__dict__}
    debug_console_hook(simulation)

# Run the simulation
t_start = time.time()
t = 0
it = 0
while t+dt <= tmax + 1e-8:
    it += 1
    t += dt
    print 'TIMESTEP %5d  -  t = %10.4f' % (it, t)
    
    update_convection(t, dt)
    
    error = 1e10
    inner_iter = 0
    while inner_iter < min_inner_iter or error > u_star_error_limit:
        inner_iter += 1
        
        momentum_prediction(t, dt)
        pressure_correction()
        velocity_update()
        pressure_update()
        
        # DEBUG:
        u0a = u_list[0].vector().array()
        u1a = u_list[1].vector().array()
        maxvel = (u0a**2+u1a**2).max()**0.5
        print '        Inner iter %3d, max velocity %15.4e' % (inner_iter, maxvel),
        error = errornorm(u_list, u_star_list)
        print ' errornorm*: %10.3e' % error,
        errorC = errornorm(u_list, u_conv_list)
        print ' errornormC: %10.3e' % errorC,
        errorG = ghia_error()
        print ' Ghia-err: %10.3e' % errorG
        
        if inner_iter == max_inner_iter:
            print '        Max inner iter reached!'
            break
    
    velocity_update_final()
    
    if it % 10 == 1:
        plotit()
    
    # Change to the second order time derivative for velocities
    time_coeffs.assign(time_coeffs2)
    
    debug_console_hook_wrapper()
    
    if maxvel > max_allowed_velocity:
        print '\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !'
        print '\n    Diverging solution!'
        print '\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n'
        break

# Print the runtime of the functions timed with the @timeit decorator
print '\n\nTIMINGS:'
tottime = time.time() - t_start
for funcname in sorted(timeit.timings):
    durations = timeit.timings[funcname]
    print '%30s total time %10.3fs for %5d runs, minimum runtime %7.3fs' % \
          (funcname, sum(durations), len(durations), min(durations)),
    print '  (%5.1f%% of tot.time)' % (sum(durations)/tottime*100) 
print '%30s total time %10.3fs\n' % ('ALL', tottime)

###########################################################################################

if PLOT:
    #dolfin.plot(marker, title='boundary')
    dolfin.plot(u_vec, title='u')
    #dolfin.plot(u_list[0], title='u_x')
    #dolfin.plot(u_list[1], title='u_y')
    #dolfin.plot(p, title='p')
    dolfin.interactive()

pyplot.ioff()
pyplot.show()
#raw_input('Press enter to quit')