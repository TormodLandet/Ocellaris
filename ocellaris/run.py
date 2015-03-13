import sys, time, collections, traceback
import yaml
import dolfin
from .multiphase import get_multi_phase_model
from .solvers import get_solver
from .boundary_conditions import BoundaryRegion
from .postprocess import setup_probes
from .utils import timeit, run_debug_console, debug_console_hook, report_error
from ocellaris.utils.code_runner import RunnablePythonString

def run_simulation(simulation):
    """
    Prepare and run a simulation
    """
    simulation.log.info('Preparing simulation ...\n')
    t_start = time.time()
    
    # Load the mesh. The mesh determines if we are in 2D or 3D
    load_mesh(simulation)
    
    # Mark the boundaries of the domain with separate marks
    # for each regions Creates a new "ds" measure
    mark_boundaries(simulation)
    
    # Load the periodic boundary conditions. This must 
    # be done before creating the function spaces as
    # they depend on the periodic constrained domain
    setup_periodic_domain(simulation)
    
    # Create function spaces. This must be done before
    # creating Dirichlet boundary conditions
    setup_function_spaces(simulation)
    
    # Load the boundary conditions. This must be done
    # before creating the solver as the solver needs
    # the Neumann conditions to define weak forms
    setup_boundary_conditions(simulation)
    
    # Setup physical constants
    g = simulation.input.get_value('physical_properties/g', required_type='list(float)')
    assert len(g) == simulation.ndim
    simulation.data['g'] = dolfin.Constant(g)
    
    # Get the density and viscosity properties from the multi phase model
    multiphase_model_name = simulation.input.get_value('multiphase_solver/type', 'SinglePhase', 'string')
    multiphase_model = get_multi_phase_model(multiphase_model_name)(simulation)
    simulation.data['rho'] = multiphase_model.get_density()
    simulation.data['nu'] = multiphase_model.get_laminar_kinematic_viscosity()
    simulation.multi_phase_model = multiphase_model
    simulation.hooks.add_pre_timestep_hook(multiphase_model.update)
    
    # Get the solver
    solver_name = simulation.input.get_value('solver/type', required_type='string')
    solver = get_solver(solver_name)(simulation)
    simulation.solver = solver
    
    # Setup postprocessing probes
    setup_probes(simulation)
    
    # Initialise the fields
    setup_initial_conditions(simulation)
    
    # Setup any hooks that may be present on the input file
    setup_hooks(simulation) 
    
    # Print information about configuration parameters
    simulation.log.info('\nPreparing simulation done in %.3f seconds' % (time.time() - t_start))
    simulation.log.info('\nSimulation configuration:')
    simulation.log.info('{:-^40}'.format(' input begin '))
    
    inp = collections.OrderedDict(simulation.input.items())
    simulation.log.info(yaml.dump(inp, indent=4))
    simulation.log.info('{:-^40}'.format(' input end '))
    simulation.log.info("\nRunning simulation ...\n")
    t_start = time.time()
    
    # Setup the debug console to optionally run at the end of each timestep
    simulation.hooks.add_post_timestep_hook(lambda report: debug_console_hook(simulation))
    
    # Setup the summary to show after the simulation
    hook = lambda success: summarise_simulation_after_running(simulation, t_start, success)
    simulation.hooks.add_post_simulation_hook(hook)
    
    # Run the simulation
    try:
        solver.run()
        success = True
    except Exception, e:
        success = False
        simulation.log.error('=== EXCEPTION =='*5)
        tb = traceback.format_tb(sys.exc_info()[2])
        simulation.log.error('Traceback:\n\n%s\n' % ''.join(tb))
        simulation.log.error('Got exception when running solver:\n%s' % str(e))
        simulation.log.error('=== EXCEPTION =='*5)
        simulation.hooks.simulation_ended(success)
    
    # Show dolfin plots?
    if simulation.input.get_value('output/plot_at_end', False, 'bool'):
        plot_at_end(simulation)
    
    # Optionally show the console for debugging and ad-hoc posprocessing
    console_at_end = simulation.input.get_value('console_at_end', False, 'bool')
    console_on_error = simulation.input.get_value('console_on_error', False, 'bool')
    if console_at_end  or (not success and console_on_error):
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


def mark_boundaries(simulation):
    """
    Mark the boundaries of the mesh with different numbers to be able to
    apply different boundary conditions to different regions 
    """
    # Create a function to mark the external facets
    marker = dolfin.FacetFunction("size_t", simulation.data['mesh'])
    
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


def setup_periodic_domain(simulation):
    """
    We need to create a constrained domain in case there are periodic 
    boundary conditions.
    """
    # This will be overwritten if there are periodic boundary conditions
    simulation.data['constrained_domain'] = None
    
    for region in simulation.data['boundary']:
        region.create_periodic_boundary_conditions()


def setup_function_spaces(simulation):
    """
    Create function spaces for velocity and pressure
    """
    # Get function space names
    Vu_name = simulation.input.get_value('solver/function_space_velocity', 'Lagrange', 'string')
    Vp_name = simulation.input.get_value('solver/function_space_pressure', 'Lagrange', 'string')
    
    # Get function space polynomial degrees
    Pu = simulation.input.get_value('solver/polynomial_degree_velocity', 1, 'int')
    Pp = simulation.input.get_value('solver/polynomial_degree_pressure', 1, 'int')
    
    # Get the constraine ddomain
    cd = simulation.data['constrained_domain']
    if cd is None:
        simulation.log.info('Creating function spaces')
    else:
        simulation.log.info('Creating function spaces with periodic boundaries')
    
    # Create the function spaces
    mesh = simulation.data['mesh']
    Vu = dolfin.FunctionSpace(mesh, Vu_name, Pu, constrained_domain=cd)
    Vp = dolfin.FunctionSpace(mesh, Vp_name, Pp, constrained_domain=cd)
    simulation.data['Vu'] = Vu
    simulation.data['Vp'] = Vp
    
    # Setup colour function if it is needed in the simulation
    multiphase_model_name = simulation.input.get_value('multiphase_solver/type', 'SinglePhase', 'string')
    if multiphase_model_name != 'SinglePhase':
        # Get the function space name and polynomial degree
        Vc_name = simulation.input.get_value('multiphase_solver/function_space_colour',
                                             'Discontinuous Lagrange', 'string')
        Pc = simulation.input.get_value('multiphase_solver/polynomial_degree_colour', 0, 'int')
        
        # Create and store function space
        Vc = dolfin.FunctionSpace(mesh, Vc_name, Pc, constrained_domain=cd)
        simulation.data['Vc'] = Vc


def setup_boundary_conditions(simulation):
    """
    Setup boundary conditions based on the simulation input
    """
    # Make dicts to gather Dirichlet and Neumann boundary conditions
    simulation.data['dirichlet_bcs'] = {}
    simulation.data['neumann_bcs'] = {}
    
    for region in simulation.data['boundary']:
        region.create_boundary_conditions()


def setup_initial_conditions(simulation):
    """
    Setup the initial values for the fields
    """
    ic = simulation.input.get_value('initial_conditions', {}, 'dict(string:dict)')
    for name, info in ic.items():
        name = str(name)
        
        if not 'p' in name:
            report_error('Invalid initial condition',
                         'You have given initial conditions for %r but this does '
                         'not seem to be a previous or pressure field.\n\n'
                         'Valid names: up0, up1, ... , p, cp, ...' % name)
        
        if not 'cpp_code' in info:
            report_error('Invalid initial condition',
                         'You have not given "cpp_code" for %r' % name)
        
        cpp_code = str(info['cpp_code'])
        
        if not name in simulation.data:
            report_error('Invalid initial condition',
                         'You have given initial conditions for %r but this does '
                         'not seem to be an existing field.' % name)
        
        func = simulation.data[name]
        V = func.function_space()
        
        available_vars = {'t': simulation.time}
        for k, v in simulation.data.items():
            if isinstance(v, (float, int, long)):
                available_vars[k] = v
            if isinstance(v, dolfin.Constant) and v.ufl_shape == ():
                available_vars[k] = v
        
        try:
            expr = dolfin.Expression(cpp_code, **available_vars)
            dolfin.project(expr, V=V, function=func)
        except Exception as e:
            report_error('Error in C++ code',
                         'The C++ code for initial conditions for %r does not compile.\n'
                         '\nCode:\n%s'
                         '\n\nError:\n%s' % (name, cpp_code, str(e)))

def setup_hooks(simulation):
    """
    Install the hooks that are given on the input file
    """
    hooks = simulation.input.get_value('hooks', {}, 'dict(string:list)')
    
    def make_hook_from_code_string(code_string, description):
        runnable = RunnablePythonString(simulation, code_string, description)
        def hook(*args, **kwargs):
            runnable.run()
        return hook
    
    hook_types = [('post_timestep', simulation.hooks.add_post_timestep_hook)]
    
    for hook_name, register_hook in hook_types:
        for hook_info in hooks.get(hook_name, []):
            name = hook_info.get('name', 'unnamed')
            enabled = hook_info.get('enabled', True)
            description = '%s hook "%s"' % (hook_name, name)
            
            if not enabled:
                simulation.log.info('Skipping disabled %s' % description)
                continue
            
            code_string = hook_info['code']
            hook = make_hook_from_code_string(code_string, description)
            register_hook(hook)
            simulation.log.info('Registering %s' % description)


def summarise_simulation_after_running(simulation, t_start, success):
    """
    Print a summary of the time spent on each part of the simulation
    """
    simulation.log.debug('\nGlobal simulation data at end of simulation:')
    for key, value in sorted(simulation.data.items()):
        simulation.log.debug('%20s = %s' % (key, repr(type(value))[:57]))
    
    # Print the runtime of the functions timed with the @timeit decorator
    simulation.log.info('\nSummary of time spent:')
    tottime = time.time() - t_start
    for funcname in sorted(timeit.timings):
        durations = timeit.timings[funcname]
        simulation.log.info('  %30s total time %7.3fs for %5d runs, minimum runtime %7.3fs'
                            % (funcname, sum(durations), len(durations), min(durations)) +
                            '  (%5.1f%% of tot.time)' % (sum(durations)/tottime*100))
    
    # Show the total duration
    h = int(tottime/60**2)
    m = int((tottime - h*60**2)/60)
    s = tottime - h*60**2 - m*60
    humantime = '%d hours %d minutes and %d seconds' % (h, m, s)
    simulation.log.info('\nSimulation done in %.3f seconds (%s)' % (tottime, humantime))


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
