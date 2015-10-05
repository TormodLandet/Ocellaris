import os
import collections
import numpy
import dolfin
from ocellaris.multiphase import get_multi_phase_model
from ocellaris import Simulation
from ocellaris.run import load_mesh, setup_boundary_conditions, mark_boundaries
from ocellaris.postprocess import setup_probes

thisdir = os.path.dirname(os.path.abspath(__file__))
inpfile = os.path.join(thisdir, 'diagonal_advection_of_square.inp')

sim = Simulation()
sim.input.read_yaml(inpfile)

VEL = numpy.array([1.0, 1.0], float)
VEL_TURN_TIME = 0.5

TMAX = 1.0
Nt = 1000
TS_MAX = 200000

PLOT = False
PNG_OUTPUT_FREQUENCY = 1
PLOT_INTERPOLATED = False
FIELDS_TO_PLOT = 'c c_beta c_grad'.split()

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

dolfin.set_log_level(dolfin.WARNING)

# Prepare simulation
mesh_facet_regions = load_mesh(sim)
mark_boundaries(sim, mesh_facet_regions)
sim.data['constrained_domain'] = None
mesh = sim.data['mesh']
sim.data['Vc'] = dolfin.FunctionSpace(mesh, "DG", 0)
setup_boundary_conditions(sim)

# Initialize the convecting velocity field
sim.data['Vu'] = dolfin.FunctionSpace(mesh, "DG", 1)
u0 = dolfin.Function(sim.data['Vu'])
u1 = dolfin.Function(sim.data['Vu'])
sim.data['u'] = sim.data['up'] = dolfin.as_vector([u0, u1])
sim.data['Vp'] = sim.data['Vc']
sim.data['p'] = dolfin.Function(sim.data['Vp'])

# Initialize the VOF scheme
multiphase_model = get_multi_phase_model('BlendedAlgebraicVOF')
vof = multiphase_model(sim)

# Initialize the colour function field
c0 = dolfin.interpolate(cexpr, sim.data['Vc'])

# Initialise probes used for reporting iso surfaces
setup_probes(sim)

################################################################################
# Runtime postprocessing

def postprocess(cfunc, t):
    cvec = cfunc.vector().array()
    
    cexpr.t = t
    target = dolfin.interpolate(cexpr, sim.data['Vc']).vector().array()
    error = ((target - cvec)**2).sum() *  60*60*2/mesh.num_cells()
    
    sim.reporting.report_timestep_value('error', error)

###############################################################################
# Time loop
print 'CFL ~', (TMAX/Nt)/mesh.hmin()

t_vec = numpy.linspace(0, TMAX, Nt)
dt_vec = t_vec[1:] - t_vec[:-1]

sim.data['cp'].assign(c0)

# Not needed, but nice to have for visualization of the 0th time step
sim.data['c'].assign(c0)
vof.convection_scheme.gradient_reconstructor.reconstruct()
    
# Dump input data
import yaml
inp = collections.OrderedDict(sim.input.items())
print yaml.dump(inp, indent=4)

postprocess(sim.data['c'], 0.0)
sim.hooks.simulation_started()
sim.reporting.log_timestep_reports()

for it in xrange(1, Nt):
    t = t_vec[it]
    dt = dt_vec[it-1]
    sim.hooks.new_timestep(it, t, dt)

    if t < VEL_TURN_TIME:
        u0.vector()[:] = VEL[0]
        u1.vector()[:] = VEL[1]
    else:
        u0.vector()[:] = -VEL[0]
        u1.vector()[:] = -VEL[1]
    
    if PLOT and it % PNG_OUTPUT_FREQUENCY == 0:
        for field in FIELDS_TO_PLOT:
            sim.plotting.plot(field)
    
    postprocess(sim.data['c'], t)
    sim.hooks.end_timestep()
    
    if it == TS_MAX:
        break

sim.hooks.simulation_ended(success=True)
sim.plotting.plot_all()
