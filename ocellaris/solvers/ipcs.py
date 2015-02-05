import dolfin
from . import Solver, register_solver

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
        
        # Register velocity functions
        uvec, upvec, uppvec = [], [], []
        for d in range(sim.ndim):
            u = dolfin.Function(Vu)
            up = dolfin.Function(Vu)
            upp = dolfin.Function(Vu)
            sim.data['u%d' % d] = u
            sim.data['up%d' % d] = up
            sim.data['upp%d' % d] = upp
            uvec.append(u)
            upvec.append(up)
            uppvec.append(upp)
        sim.data['u'] = uvec = dolfin.as_vector(uvec)
        sim.data['up'] = upvec = dolfin.as_vector(upvec)
        sim.data['upp'] = uppvec = dolfin.as_vector(uppvec)
        
        # Register pressure functions
        sim.data['p'] = p = dolfin.Function(Vp)

    def velocity_prediction(self):
        pass

    def pressure_correction(self):
        pass
    
    def run(self):
        """
        Run the simulation
        """
        dt = self.simulation.input['time']['dt']
        tmax = self.simulation.input['time']['tmax']
