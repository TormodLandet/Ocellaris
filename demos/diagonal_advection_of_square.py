import time
import numpy
import dolfin
from dgvof.vof import BlendedAlgebraicVofScheme
from dgvof.plot import Plot2DDG0, plot_2d_DG0
from tictoc import tic, toc

dolfin.set_log_level(dolfin.WARNING)

xmax = 1.5; Nx = 60
ymax = 1.5; Ny = 60
VEL = numpy.array([1.0, 1.0], float)
VEL_TURN_TIME = 0.5
TMAX = 1.0
Nt = 300
HRIC_FORCE_UPWIND = False
PLOT = False
PLOT_INTERPOLATED = True
PNG_OUTPUT_FREQUENCY = 1

print 'CFL ~', (TMAX/Nt)/(xmax/Nx)*VEL[0], (TMAX/Nt)/(ymax/Ny)*VEL[1]

class InitialC(dolfin.Expression):
    def eval(self, value, x):
        value[0] = 1 if (0.25 < x[0] < 0.50 and 0.25 < x[1] < 0.50) else 0
c0expr = InitialC()

################################################################################
# Problem definition

# Initialize mesh and function space
mesh = dolfin.RectangleMesh(0, 0, xmax, ymax, Nx, Ny, 'left/right')

# Initialize the convecting velocity field
vel_func_space = dolfin.VectorFunctionSpace(mesh, "DG", 0)
vel = dolfin.Function(vel_func_space)

# Initialize the VOF scheme
vof = BlendedAlgebraicVofScheme(mesh, 'HRIC', vel)
vof.convection_scheme.force_upwind = HRIC_FORCE_UPWIND

# Initialize the colour function field
c0 = dolfin.interpolate(c0expr, vof.function_space)

################################################################################
# Runtime postprocessing

class RuntimeOutput(object):
    def __init__(self):
        self.prevtime = time.time()
        
    def __call__(self, t, cfunc):
        cvec = cfunc.vector().array()
        
        csum = dolfin.assemble(cfunc*dolfin.dx)
        cmax = cvec.max()
        cmin = cvec.min()
    
        now = time.time()
        runtime = now - self.prevtime
        print "Time = %6.3f - runtime = %5.3f - csum = %8.5f - minmax = %8.5f %8.5f" % (t, runtime, csum, cmin, cmax)
        self.prevtime = now

###############################################################################
# Time loop

t_vec = numpy.linspace(0, TMAX, Nt)
dt_vec = t_vec[1:] - t_vec[:-1]

def make_movie_frame():
    movie_png_filename = "fig/cfunc_%05d_%010.5f.png" % (it, t)
    print movie_png_filename
    movie.plot(movie_png_filename, skip_zero_values=False)
    
    plot_2d_DG0(movie.verts, vof.convection_scheme.gradient_reconstructor.gradient[:,0],
                "fig/gradient_x_cfunc_%05d_%010.5f.png" % (it, t),
                xlim=(movie.coords[:,0].min(), movie.coords[:,0].max()),
                ylim=(movie.coords[:,1].min(), movie.coords[:,1].max()),
                cmap='seismic')
    
    plot_2d_DG0(movie.verts, vof.convection_scheme.gradient_reconstructor.gradient[:,1],
                "fig/gradient_y_cfunc_%05d_%010.5f.png" % (it, t),
                xlim=(movie.coords[:,0].min(), movie.coords[:,0].max()),
                ylim=(movie.coords[:,1].min(), movie.coords[:,1].max()),
                cmap='seismic')

vof.prev_colour_function.assign(c0)

# Not needed, but nice to have for visualization of the 0th time step
vof.colour_function.assign(c0)
vof.convection_scheme.gradient_reconstructor.reconstruct()

# Function object dumping interesting variables to screen
runtime_output = RuntimeOutput()
runtime_output(t_vec[0], vof.colour_function)

# Make png frames of the evolution of the colour function
movie = Plot2DDG0(vof.colour_function)
t = t_vec[0]; it = 0
make_movie_frame()

tic('timeloop')
for it in xrange(1, Nt):
    t = t_vec[it]
    dt = dt_vec[it-1]

    if t < VEL_TURN_TIME:
        vel.assign(dolfin.Constant(VEL))
    else:
        vel.assign(dolfin.Constant(-VEL))

    vof.update(t, dt)
    runtime_output(t, vof.colour_function)
    if it % PNG_OUTPUT_FREQUENCY == 0:
        make_movie_frame()
    
    #plot(c_face, title='c_f at t=%f'%t, wireframe=True)
    if PLOT and it % 30 == 0 and it != 0:
        dolfin.plot(vof.colour_function, title='c at t=%f'%t)
    #interactive()
    #if it > 2: break
toc()

###############################################################################
# Visualisation

if PLOT_INTERPOLATED:
    # Create a CG mesh and function space for visualisation
    mesh_vis = dolfin.RectangleMesh(0, 0, xmax, ymax, Nx*5, Ny*5, 'left/right')
    V_vis = dolfin.FunctionSpace(mesh_vis, "CG", 1)

    # Visualize the initial colour field
    c0_vis = dolfin.project(c0, V=V_vis)
    dolfin.plot(c0_vis, title='c - initial', wireframe=False)

    # Visualize the final colour field
    c_vis = dolfin.project(vof.colour_function, V=V_vis)
    dolfin.plot(c_vis, title='c - final', wireframe=False)

    #plot(mesh, wireframe=True, title='Mesh DG')

if PLOT or PLOT_INTERPOLATED:
    dolfin.interactive()
