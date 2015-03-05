import time
import json
import dolfin
from .multiphase import get_multi_phase_model
from .solvers import get_solver
from .boundary_conditions import BoundaryRegion
from .postprocess import setup_probes
from .utils import timeit, run_debug_console, debug_console_hook

def run_simulation(simulation):
    """
    Prepare and run a simulation
    """
    # Set log levels
    available_log_levels = {'critical': dolfin.CRITICAL,
                            'error': dolfin.ERROR,
                            'warning': dolfin.WARNING,
                            'info': dolfin.INFO,
                            'progress': dolfin.PROGRESS,
                            'debug': dolfin.DEBUG}
    # Ocellaris log level
    log_level = simulation.input.get_value('output/ocellaris_log_level', 'info')
    simulation.log.set_log_level(available_log_levels[log_level])
    # Dolfin log level
    df_log_level = simulation.input.get_value('output/dolfin_log_level', 'warning')
    dolfin.set_log_level(available_log_levels[df_log_level])
    
    simulation.log.info('Preparing simulation ...\n')
    t1 = time.time()
    
    # Load the mesh
    load_mesh(simulation)
    
    # Create function spaces. This must be done
    # early as it is needed by the DirichletBC
    make_function_spaces(simulation)
    
    # Load the boundary conditions. This must be done
    # before creating the solver as the solver needs
    # the Neumann conditions to define weak forms
    setup_boundary_conditions(simulation)
    
    # Setup physical constants
    g = simulation.input.get_value('physical_properties/g', required_type='list(float)')
    assert len(g) == simulation.ndim
    simulation.data['g'] = dolfin.Constant(g)
    
    # Get the density and viscosity properties from the multi phase model
    multiphase_model_name = simulation.input.get_value('multiphase_model', 'SinglePhase', 'string')
    multiphase_model = get_multi_phase_model(multiphase_model_name)(simulation)
    simulation.data['rho'] = multiphase_model.get_density()
    simulation.data['nu'] = multiphase_model.get_laminar_kinematic_viscosity()
    simulation.add_pre_timestep_hook(multiphase_model.update)
    
    # Get the solver
    solver_name = simulation.input.get_value('solver/type', required_type='string')
    solver = get_solver(solver_name)(simulation)
    simulation.data['solver'] = solver
    
    # Setup postprocessing probes
    setup_probes(simulation)
    
    # Print information about configuration parameters
    simulation.log.info('\nPreparing simulation done in %.3f seconds' % (time.time() - t1))
    simulation.log.info('\nSimulation configuration:')
    simulation.log.info('{:-^40}'.format(' input begin '))
    simulation.log.info(json.dumps(simulation.input, indent=4))
    simulation.log.info('{:-^40}'.format(' input end '))
    simulation.log.info("\nRunning simulation ...\n")
    t1 = time.time()
    
    # Setup the debug console to optionally run at the end of each timestep
    simulation.add_post_timestep_hook(lambda report: debug_console_hook(simulation))
    
    # Run the simulation
    solver.run()
    
    simulation.log.debug('\nGlobal simulation data at end of simulation:')
    for key, value in sorted(simulation.data.items()):
        simulation.log.debug('%20s = %s' % (key, repr(type(value))[:57]))
    
    # Print the runtime of the functions timed with the @timeit decorator
    simulation.log.info('\nSummary of time spent:')
    tottime = time.time() - t1
    for funcname in sorted(timeit.timings):
        durations = timeit.timings[funcname]
        simulation.log.info('  %30s total time %7.3fs for %5d runs, minimum runtime %7.3fs'
                            % (funcname, sum(durations), len(durations), min(durations)) +
                            '  (%5.1f%% of tot.time)' % (sum(durations)/tottime*100))
    simulation.log.info('\nSimulation done in %.3f seconds' % tottime)
    
    if simulation.input.get_value('output/plot_at_end', False, 'bool'):
        plot_at_end(simulation)
    
    if simulation.input.get_value('console_at_end', False, 'bool'):
        run_debug_console(simulation)

def load_mesh(simulation):
    """
    Get the mesh from the simulation input
    """
    inp = simulation.input
    mesh_type = inp['mesh']['type']
    
    if mesh_type == 'Rectangle':
        startx = inp['mesh'].get('startx', 0)
        starty = inp['mesh'].get('starty', 0)
        endx = inp['mesh'].get('endx', 1)
        endy = inp['mesh'].get('endy', 1)
        Nx = int(inp['mesh']['Nx'])
        Ny = int(inp['mesh']['Ny'])
        diagonal = inp['mesh'].get('diagonal', 'left/right')
        mesh = dolfin.RectangleMesh(startx, starty, endx, endy, Nx, Ny, diagonal)
    
    simulation.set_mesh(mesh)

def make_function_spaces(simulation):
    """
    Create function spaces for velocity and pressure
    """
    # Get function space names
    Vu_name = simulation.input.get_value('solver/function_space_velocity', 'Lagrange', 'string')
    Vp_name = simulation.input.get_value('solver/function_space_pressure', 'Lagrange', 'string')
    
    # Make sure strings are not Python 2 unicode objects
    Vu_name = str(Vu_name)
    Vp_name = str(Vp_name)
    
    # Get function space polynomial degrees
    Pu = simulation.input.get_value('solver/polynomial_degree_velocity', 1, 'int')
    Pp = simulation.input.get_value('solver/polynomial_degree_pressure', 1, 'int')
    
    # Create the function spaces
    mesh = simulation.data['mesh']
    Vu = dolfin.FunctionSpace(mesh, Vu_name, Pu)
    Vp = dolfin.FunctionSpace(mesh, Vp_name, Pp)
    simulation.data['Vu'] = Vu
    simulation.data['Vp'] = Vp

def setup_boundary_conditions(simulation):
    """
    Setup boundary conditions based on the simulation input
    """
    # Create a function to mark the external facets
    marker = dolfin.FacetFunction("size_t", simulation.data['mesh'])
    
    # Make dicts to gather Dirichlet and Neumann boundary conditions
    simulation.data['dirichlet_bcs'] = {}
    simulation.data['neumann_bcs'] = {}
     
    # Create boundary regions and let them mark the part of the
    # boundary that they belong to. They also create boundary
    # condition objects that are later used in the eq. solvers
    boundary = []
    for index, _ in enumerate(simulation.input.get_value('boundary_conditions', required_type='list(dict)')):
        part = BoundaryRegion(simulation, marker, index)
        boundary.append(part)
    simulation.data['boundary'] = boundary
    simulation.data['boundary_marker'] = marker
    
    # Create a boundary measure that is aware of the marked regions
    ds = dolfin.Measure('ds')[marker]
    simulation.data['ds'] = ds

def plot_at_end(simulation):
    """
    Open dolfin plotting windows with results at the end of
    the simulation
    """
    # Plot velocity components
    for d in range(simulation.ndim):
        name = 'u%d' % d
        dolfin.plot(simulation.data[name], title=name)
        #name = 'u_star%d' % d
        #dolfin.plot(simulation.data[name], title=name)
    
    dolfin.plot(simulation.data['u'], title='u')
    
    # Plot pressure
    dolfin.plot(simulation.data['p'], title='p')
    
    # Plot colour function
    if 'c' in simulation.data:
        dolfin.plot(simulation.data['c'], title='c')
        
    dolfin.plot(simulation.data['boundary_marker'], title='boundary_marker')
    
    dolfin.interactive()
