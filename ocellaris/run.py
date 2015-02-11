import time
import json
import dolfin
from .multiphase import get_multi_phase_model
from .solvers import get_solver
from .boundary_conditions import BoundaryRegion

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
    log_level = simulation.input.get('output', {}).get('ocellaris_log_level', 'info')
    simulation.log.set_log_level(available_log_levels[log_level])
    # Dolfin log level
    df_log_level = simulation.input.get('output', {}).get('dolfin_log_level', 'warning')
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
    g = simulation.input['physical_properties']['g']
    assert len(g) == simulation.ndim
    simulation.data['g'] = dolfin.Constant(g)
    
    # Get the density and viscosity properties from the multi phase model
    multiphase_model_name = simulation.input.get('multiphase_model', 'SinglePhase')
    multiphase_model = get_multi_phase_model(multiphase_model_name)(simulation)
    simulation.data['rho'] = multiphase_model.get_density()
    simulation.data['nu'] = multiphase_model.get_laminar_kinematic_viscosity()
    simulation.add_pre_timestep_hook(multiphase_model.update)
    
    # Get the solver
    solver_name = simulation.input['solver']['type']
    solver = get_solver(solver_name)(simulation)
    
    # Print information about configuration parameters
    simulation.log.info('\nPreparing simulation done in %.3f seconds' % (time.time() - t1))
    simulation.log.info('\nSimulation configuration:')
    simulation.log.info('{:-^40}'.format(' input begin '))
    simulation.log.info(json.dumps(simulation.input, indent=4))
    simulation.log.info('{:-^40}'.format(' input end '))
    simulation.log.info("\nRunning simulation ...\n")
    t1 = time.time()

    # Run the simulation 
    solver.run()
    
    simulation.log.debug('\nGlobal simulation data at end of simulation:')
    for key, value in sorted(simulation.data.items()):
        simulation.log.debug('%20s = %s' % (key, repr(type(value))[:57]))
    
    simulation.log.info('\nSimulation done in %.3f seconds' % (time.time() - t1))
        
    if simulation.input.get('output', {}).get('plot_at_end', False):
        plot_at_end(simulation)

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
    Pu = simulation.input['solver'].get('polynomial_degree_velocity', 1)
    Pp = simulation.input['solver'].get('polynomial_degree_pressure', 1)
    mesh = simulation.data['mesh']
    Vu = dolfin.FunctionSpace(mesh, 'Discontinuous Lagrange', Pu)
    Vp = dolfin.FunctionSpace(mesh, 'Discontinuous Lagrange', Pp)
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
    for index, _ in enumerate(simulation.input['boundary_conditions']):
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
    
    # Plot pressure
    dolfin.plot(simulation.data['p'], title='p')
    
    # Plot colour function
    if 'c' in simulation.data:
        dolfin.plot(simulation.data['c'], title='c')
        
    dolfin.plot(simulation.data['boundary_marker'], title='boundary_marker')
    
    dolfin.interactive()
