import time
import numpy
import dolfin
from dgvof.vof import BlendedAlgebraicVofScheme
from dgvof import Simulation

dolfin.set_log_level(dolfin.WARNING)

xmax = 1.5; Nx = 128
ymax = 1.5; Ny = 128

VEL = numpy.array([1.0, 1.0], float)
VEL_TURN_TIME = 0.5

TMAX = 1.0
Nt = 1000
TS_MAX = 20

PLOT = False
PNG_OUTPUT_FREQUENCY = 10
PLOT_INTERPOLATED = False


print 'CFL ~', (TMAX/Nt)/(xmax/Nx)*VEL[0], (TMAX/Nt)/(ymax/Ny)*VEL[1]

class CField(dolfin.Expression):
    t = 0
    def eval(self, value, x):
        if self.t < VEL_TURN_TIME:
            cx = 0.375 + VEL[0]*self.t
            cy = 0.375 + VEL[1]*self.t
        else:
            cx = 0.375 + 2*VEL[0]*VEL_TURN_TIME - VEL[0]*self.t
            cy = 0.375 + 2*VEL[1]*VEL_TURN_TIME - VEL[1]*self.t
        value[0] = 1 if (cx-0.125 < x[0] < cx+0.125 and cy-0.125 < x[1] < cy+0.125) else 0
cexpr = CField()

################################################################################
# Problem definition

sim = Simulation()
sim.read_json_input_file('input.json')

# Initialize mesh
mesh = dolfin.RectangleMesh(0, 0, xmax, ymax, Nx, Ny, 'left/right')
sim.set_mesh(mesh)

# Initialize the convecting velocity field
vel_func_space = dolfin.VectorFunctionSpace(mesh, "DG", 1)
vel = dolfin.Function(vel_func_space)
sim.data['u'] = vel

# Initialize the VOF scheme
vof = BlendedAlgebraicVofScheme(sim)

# Initialize the colour function field
c0 = dolfin.interpolate(cexpr, vof.function_space)

################################################################################
# Runtime postprocessing

def postprocess(cfunc, t):
    cvec = cfunc.vector().array()
    
    csum = dolfin.assemble(cfunc*dolfin.dx)
    cmax = cvec.max()
    cmin = cvec.min()
    
    cexpr.t = t
    target = dolfin.interpolate(cexpr, vof.function_space).vector().array()
    error = ((target - cvec)**2).sum() *  60*60/(Nx*Ny)
    
    sim.reporting.report_timestep_value('csum', csum)
    sim.reporting.report_timestep_value('cmin', cmin)
    sim.reporting.report_timestep_value('cmax', cmax)
    sim.reporting.report_timestep_value('error', error)
        

###############################################################################
# Time loop

t_vec = numpy.linspace(0, TMAX, Nt)
dt_vec = t_vec[1:] - t_vec[:-1]

vof.prev_colour_function.assign(c0)

# Not needed, but nice to have for visualization of the 0th time step
vof.colour_function.assign(c0)
vof.convection_scheme.gradient_reconstructor.reconstruct()


# Dump input data
import json
print json.dumps(sim.input, indent=4)

# Make png frames of the evolution of the colour function
sim.plotting.plot_all()

postprocess(vof.colour_function, 0.0)
sim.reporting.print_timestep_report()

for it in xrange(1, Nt):
    t = t_vec[it]
    dt = dt_vec[it-1]
    sim.new_timestep(it, t, dt)

    if t < VEL_TURN_TIME:
        vel.assign(dolfin.Constant(VEL))
    else:
        vel.assign(dolfin.Constant(-VEL))

    vof.update(t, dt)
    
    if PLOT and it % PNG_OUTPUT_FREQUENCY == 0:
        sim.plotting.plot_all()
    
    postprocess(vof.colour_function, t)
    sim.end_timestep()
    
    if it == TS_MAX:
        break

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

if PLOT_INTERPOLATED:
    dolfin.interactive()

sim.plotting.plot_all()