import sys, time, importlib
import dolfin
from ocellaris.solvers import get_solver
from ocellaris.probes import setup_probes
from ocellaris.utils import interactive_console_hook, ocellaris_error, ocellaris_interpolate, log_timings
from ocellaris.utils import RunnablePythonString, OcellarisCppExpression
from ocellaris.utils import verify_key
from ocellaris.solver_parts import BoundaryRegion, get_multi_phase_model, MeshMorpher


def setup_simulation(simulation):
    """
    Prepare an Ocellaris simulation for running
    """
    # Setup user code (custom sys.path, custom module imports)
    setup_user_code(simulation)
    
    simulation.log.info('Preparing simulation')
    simulation.log.info('Output prefix is: %s' % 
                        simulation.input.get_value('output/prefix', '', 'string'))
    simulation.log.info("Current time: %s\n" % time.strftime('%Y-%m-%d %H:%M:%S'))
    t_start = time.time()
    
    # Set linear algebra backend to PETSc
    setup_fenics(simulation)
    
    # Make time and timestep available in expressions for the initial conditions etc
    simulation.log.info('Creating time simulation')
    simulation.time = simulation.input.get_value('time/tstart', 0.0, 'float')
    simulation.dt = simulation.input.get_value('time/dt', required_type='float')
    assert simulation.dt > 0
    
    # Preliminaries, before setup begins
    simulation.flush()
    
    # Get the multi phase model class
    multiphase_model_name = simulation.input.get_value('multiphase_solver/type', 'SinglePhase', 'string')
    multiphase_class = get_multi_phase_model(multiphase_model_name)
    
    # Get the main solver class
    solver_name = simulation.input.get_value('solver/type', required_type='string')
    simulation.log.info('Creating equation solver %s' % solver_name)
    solver_class = get_solver(solver_name)
    
    ###########################################################################
    # Setup the Ocellaris simulation
    
    if not simulation.restarted:
        # Load the mesh. The mesh determines if we are in 2D or 3D
        load_mesh(simulation)
    
    # Mark the boundaries of the domain with separate marks
    # for each regions Creates a new "ds" measure
    mark_boundaries(simulation)
    simulation.flush()
    
    # Load the periodic boundary conditions. This must 
    # be done before creating the function spaces as
    # they depend on the periodic constrained domain
    setup_periodic_domain(simulation)
    simulation.flush()
    
    # Create function spaces. This must be done before
    # creating Dirichlet boundary conditions
    setup_function_spaces(simulation, solver_class, multiphase_class)
    simulation.flush()
    
    # Load the mesh morpher used for prescribed mesh velocities and ALE multiphase solvers
    simulation.mesh_morpher = MeshMorpher(simulation)
    simulation.flush()
    
    # Setup physical constants and multi-phase model (g, rho, nu, mu)
    setup_physical_properties(simulation, multiphase_class)
    simulation.flush()
    
    # Load the boundary conditions. This must be done
    # before creating the solver as the solver needs
    # the Neumann conditions to define weak forms
    setup_boundary_conditions(simulation)
    simulation.flush()
            
    # Create momentum sources (usefull for MMS tests etc)
    setup_sources(simulation)
    simulation.flush()
    
    # Create the solver
    simulation.solver = solver_class(simulation)
    
    # Setup postprocessing probes
    setup_probes(simulation)
    simulation.flush()
    
    # Initialise the fields
    if not simulation.restarted:
        setup_initial_conditions(simulation)
        simulation.flush()
        
    # Setup the solution properties
    simulation.solution_properties.setup()
    
    # Setup any hooks that may be present on the input file
    setup_hooks(simulation)
    
    # Setup the interactive console to optionally run at the end of each timestep
    simulation.hooks.add_post_timestep_hook(lambda: interactive_console_hook(simulation),
                                            'Interactive console commands')
    
    # Setup the summary to show after the simulation
    hook = lambda success: summarise_simulation_after_running(simulation, success)
    simulation.hooks.add_post_simulation_hook(hook, 'Summarise simulation')
    
    # Show time spent setting up the solver
    simulation.log.info('\nPreparing simulation done in %.3f seconds' % (time.time() - t_start))
    
    # Show all registered hooks
    simulation.hooks.show_hook_info()
    simulation.flush()


def setup_user_code(simulation):
    """
    Setup custom path and user code imports
    """
    # Extend the Python search path
    paths = simulation.input.get_value('user_code/python_path', [], 'list(string)')
    for path in paths[::-1]:
        sys.path.insert(0, path)
    
    # Import modules, possibly defining new solvers etc
    modules = simulation.input.get_value('user_code/modules', [], 'list(string)')
    for module in modules:
        importlib.import_module(module.strip())
        
    # Code to run right before setup (custom mesh generation etc)
    code = simulation.input.get_value('user_code/code', '', 'string')
    if code:
        consts = simulation.input.get_value('user_code/constants', {}, 'dict(string:any)')
        variables = consts.copy()
        variables['simulation'] = simulation
        variables['__file__'] = simulation.input.file_name
        try:
            exec(code, globals(), variables)
        except Exception:
            simulation.log.error('Exception in input file user_code/code')
            raise


def setup_fenics(simulation):
    """
    Setup FEniCS parameters like linear algebra backend
    """
    # Test for PETSc linear algebra backend
    if not dolfin.has_linear_algebra_backend("PETSc"):
        ocellaris_error('Missing PETSc',
                        'DOLFIN has not been configured with PETSc '
                        'which is needed by Ocellaris.')
    dolfin.parameters['linear_algebra_backend'] = 'PETSc'
    
    # Form compiler "uflacs" needed for isoparametric elements
    form_compiler = simulation.input.get_value('solver/form_compiler', 'auto', 'string')
    dolfin.parameters['form_compiler']['representation'] = form_compiler


def load_mesh(simulation):
    """
    Get the mesh from the simulation input
    
    Returns the facet regions contained in the mesh data
    or None if these do not exist
    """
    inp = simulation.input
    mesh_type = inp.get_value('mesh/type', required_type='string')
    mesh_facet_regions = None
    
    # For testing COMM_SELF may be specified
    comm_type = inp.get_value('mesh/mpi_comm', 'WORLD', required_type='string')
    verify_key('mesh/mpi_comm', comm_type, ('SELF', 'WORLD'))
    if comm_type == 'WORLD':
        comm = dolfin.MPI.comm_world
    else:
        comm = dolfin.MPI.comm_self
    
    if mesh_type == 'Rectangle':
        simulation.log.info('Creating rectangular mesh')
        
        startx = inp.get_value('mesh/startx', 0, 'float')
        starty = inp.get_value('mesh/starty', 0, 'float')
        start = dolfin.Point(startx, starty)
        endx = inp.get_value('mesh/endx', 1, 'float')
        endy = inp.get_value('mesh/endy', 1, 'float')
        end = dolfin.Point(endx, endy)
        Nx = inp.get_value('mesh/Nx', required_type='int')
        Ny = inp.get_value('mesh/Ny', required_type='int')
        diagonal = inp.get_value('mesh/diagonal', 'right', required_type='string')
        
        mesh = dolfin.RectangleMesh(comm, start, end, Nx, Ny, diagonal)
    
    elif mesh_type == 'Box':
        simulation.log.info('Creating box mesh')
        
        startx = inp.get_value('mesh/startx', 0, 'float')
        starty = inp.get_value('mesh/starty', 0, 'float')
        startz = inp.get_value('mesh/startz', 0, 'float')
        start = dolfin.Point(startx, starty, startz)
        endx = inp.get_value('mesh/endx', 1, 'float')
        endy = inp.get_value('mesh/endy', 1, 'float')
        endz = inp.get_value('mesh/endz', 1, 'float')
        end = dolfin.Point(endx, endy, endz)
        Nx = inp.get_value('mesh/Nx', required_type='int')
        Ny = inp.get_value('mesh/Ny', required_type='int')
        Nz = inp.get_value('mesh/Nz', required_type='int')
        
        mesh = dolfin.BoxMesh(comm, start, end, Nx, Ny, Nz)
    
    elif mesh_type == 'UnitDisc':
        simulation.log.info('Creating circular mesh')
        
        N = inp.get_value('mesh/N', required_type='int')
        degree = inp.get_value('mesh/degree', 1, required_type='int')
        gdim = inp.get_value('mesh/gdim', 2, required_type='int')
        
        mesh = dolfin.UnitDiscMesh(comm, N, degree, gdim)
        
        if degree > 1 and dolfin.parameters['form_compiler']['representation'] != 'uflacs':
            simulation.log.warning('Using isoparametric elements without uflacs!')
    
    elif mesh_type == 'XML':
        simulation.log.info('Creating mesh from XML file')
        
        mesh_file = inp.get_value('mesh/mesh_file', required_type='string')
        facet_region_file = inp.get_value('mesh/facet_region_file', None, required_type='string')
        
        # Load the mesh from file
        pth = inp.get_input_file_path(mesh_file)
        mesh = dolfin.Mesh(comm, pth)
        
        # Load the facet regions if available
        if facet_region_file is not None:
            pth = inp.get_input_file_path(facet_region_file)
            mesh_facet_regions = dolfin.MeshFunction('size_t', mesh, pth)
        else:
            mesh_facet_regions = None
    
    elif mesh_type == 'HDF5':
        simulation.log.info('Creating mesh from DOLFIN HDF5 file')
        h5_file_name = inp.get_value('mesh/mesh_file', required_type='string')
        
        with dolfin.HDF5File(comm, h5_file_name, 'r') as h5:
            # Read mesh
            mesh = dolfin.Mesh()
            h5.read(mesh, '/mesh', False)
            
            # Read facet regions
            if h5.has_dataset('/mesh_facet_regions'):
                mesh_facet_regions = dolfin.FacetFunction('size_t', mesh)
                h5.read(mesh_facet_regions, '/mesh_facet_regions')
            else:
                mesh_facet_regions = None
    
    simulation.set_mesh(mesh, mesh_facet_regions)


def mark_boundaries(simulation):
    """
    Mark the boundaries of the mesh with different numbers to be able to
    apply different boundary conditions to different regions 
    """
    simulation.log.info('Creating boundary regions')
    
    # Create a function to mark the external facets
    mesh = simulation.data['mesh']
    marker = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    mesh_facet_regions = simulation.data['mesh_facet_regions']
    
    # Create boundary regions and let them mark the part of the
    # boundary that they belong to. They also create boundary
    # condition objects that are later used in the eq. solvers
    boundary = []
    for index, _ in enumerate(simulation.input.get_value('boundary_conditions', required_type='list(dict)')):
        part = BoundaryRegion(simulation, marker, index, mesh_facet_regions)
        boundary.append(part)
    
    simulation.data['boundary'] = boundary
    simulation.data['boundary_marker'] = marker
    simulation.data['boundary_by_name'] = {b.name: b for b in boundary}
    
    # Create a boundary measure that is aware of the marked regions
    mesh = simulation.data['mesh']
    ds = dolfin.Measure('ds', domain=mesh, subdomain_data=marker)
    simulation.data['ds'] = ds
    
    # Show region sizes
    one = dolfin.Constant(1)
    for region in boundary:
        length = dolfin.assemble(one*ds(region.mark_id, domain=mesh))
        pf = simulation.log.info if length > 0.0 else simulation.log.warning
        pf('    Boundary region %s has size %f' % (region.name, length))
    length0 = dolfin.assemble(one*ds(0, domain=mesh))
    pf = simulation.log.info if length0 == 0.0 else simulation.log.warning
    pf('    Boundary region UNMARKED has size %f' % length0)


def setup_periodic_domain(simulation):
    """
    We need to create a constrained domain in case there are periodic 
    boundary conditions.
    """
    simulation.log.info('Creating periodic boundary conditions (if specified)')
    
    # This will be overwritten if there are periodic boundary conditions
    simulation.data['constrained_domain'] = None
    
    for region in simulation.data['boundary']:
        region.create_periodic_boundary_conditions()


def setup_function_spaces(simulation, solver_class, multiphase_class):
    """
    Create function spaces for equation solver and multiphase solver
    """
    # Create pressure, velocity etc spaces
    solver_class.create_function_spaces(simulation)
    
    # Create multi phase related function space (typically density or colour function)
    multiphase_class.create_function_space(simulation)
    
    for name, V in simulation.data.items():
        if isinstance(V, dolfin.FunctionSpace):
            family = V.ufl_element().family()
            degree = V.ufl_element().degree()
            simulation.log.info('    Function space %s has dimension %d (%s degree %d)' % 
                                (name, V.dim(), family, degree))


def setup_physical_properties(simulation, multiphase_class):
    """
    Gravity vector and rho/nu/mu fields are created here
    """
    ndim = simulation.ndim
    g = simulation.input.get_value('physical_properties/g', [0]*ndim, required_type='list(float)')
    assert len(g) == simulation.ndim
    simulation.data['g'] = dolfin.Constant(g)
    
    # Get the density and viscosity properties from the multi phase model
    simulation.multi_phase_model = multiphase_class(simulation)
    
    simulation.data['rho'] = simulation.multi_phase_model.get_density(0)
    simulation.data['nu'] = simulation.multi_phase_model.get_laminar_kinematic_viscosity(0)
    simulation.data['mu'] = simulation.multi_phase_model.get_laminar_dynamic_viscosity(0)


def setup_boundary_conditions(simulation):
    """
    Setup boundary conditions based on the simulation input
    """
    simulation.log.info('Creating boundary conditions')
    
    # Make dicts to gather Dirichlet and Neumann boundary conditions
    simulation.data['dirichlet_bcs'] = {}
    simulation.data['neumann_bcs'] = {}
    simulation.data['robin_bcs'] = {}
    simulation.data['slip_bcs'] = {}
    simulation.data['outlet_bcs'] = []
    
    for region in simulation.data['boundary']:
        region.create_boundary_conditions()


def setup_sources(simulation):
    """
    Setup the momentum sources
    """
    simulation.log.info('Creating sources')
    
    ms = simulation.input.get_value('momentum_sources', {}, 'list(dict)')
    sources = []
    for i in range(len(ms)):
        name = simulation.input.get_value('momentum_sources/%d/name' % i, required_type='string')
        degree = simulation.input.get_value('momentum_sources/%d/degree' % i, required_type='int')
        cpp_code = simulation.input.get_value('momentum_sources/%d/cpp_code' % i, required_type='list(string)')
        description = 'momentum source %r' % name
        simulation.log.info('    C++ %s' % description)
        
        expr = OcellarisCppExpression(simulation, cpp_code, description, degree, update=True)
        sources.append(expr)
    simulation.data['momentum_sources'] = sources


def setup_initial_conditions(simulation):
    """
    Setup the initial values for the fields
    
    NOTE: this is never run on simulation restarts!
    """
    simulation.log.info('Creating initial conditions')
    
    ic = simulation.input.get_value('initial_conditions', {}, 'dict(string:dict)')
    has_file = False
    for name, info in ic.items():
        name = str(name)
        
        if name == 'file':
            has_file = True
            continue
        
        if not 'p' in name:
            ocellaris_error('Invalid initial condition',
                            'You have given initial conditions for %r but this does '
                            'not seem to be a previous or pressure field.\n\n'
                            'Valid names: up0, up1, ... , p, cp, rho_p, ...' % name)
        
        if not 'cpp_code' in info:
            ocellaris_error('Invalid initial condition',
                            'You have not given "cpp_code" for %r' % name)
        
        cpp_code = str(info['cpp_code'])
        
        if not name in simulation.data:
            ocellaris_error('Invalid initial condition',
                            'You have given initial conditions for %r but this does '
                            'not seem to be an existing field.' % name)
        
        func = simulation.data[name]
        V = func.function_space()
        description = 'initial conditions for %r' % name
        simulation.log.info('    C++ %s' % description)
        
        # Run the C++ code to set the initial value of the function
        ocellaris_interpolate(simulation, cpp_code, description, V, func)
    
    # Some fields start out as copies, we do that here so that the input file
    # does not have to contain superfluous initial conditions
    comp_name_pairs = [('up_conv%d', 'up%d'), ('upp_conv%d', 'upp%d')]
    for cname_pattern, cname_main_pattern in comp_name_pairs:
        for d in range(simulation.ndim):
            cname = cname_pattern % d
            cname_main = cname_main_pattern % d
            
            if cname in ic:
                simulation.log.info('    Leaving %s as set by initial condition' % cname)
                continue
            
            if cname not in simulation.data or cname_main not in ic:
                continue
            
            simulation.data[cname].assign(simulation.data[cname_main])
            simulation.log.info('    Assigning initial value %s = %s' % (cname, cname_main))
    
    if has_file:
        setup_initial_conditions_from_restart_file(simulation)


def setup_initial_conditions_from_restart_file(simulation):
    """
    Read initial function values from a restart file
    
    This code is used when starting from an input file using
    initial conditions (but nothing else) from a restart file
    
    NOTE: this is never run on simulation restarts!
    """
    inp = simulation.input.get_value('initial_conditions/file', required_type='Input')
    h5_file_name = inp.get_value('h5_file', required_type='string')
    same_mesh = inp.get_value('same_mesh', False, required_type='bool')
    simulation.log.info('    Loading data from h5 file %r' % h5_file_name)
    simulation.log.info('    Forcing same mesh: %r' % same_mesh)
    funcs = simulation.io.load_restart_file_functions(h5_file_name)
    simulation.log.info('    Found %d functions' % len(funcs))
    
    tmp = [f for f in funcs.values() if f is not None]
    if not tmp:
        simulation.log.warning('    Found no supported functions!')
        return
    
    for name, f0 in sorted(funcs.items()):
        if f0 is None:
            simulation.log.warning('    Skipping unsupported function %s' % name)
            continue
        elif name not in simulation.data:
            simulation.log.warning('    Found initial condition for %r in h5 file, but '
                                   'this function is not present in the simulation' % name)
            continue
        
        f = simulation.data[name]
        if same_mesh:
            arr = f0.vector().get_local()
            f.vector().set_local(arr)
            f.vector().apply('insert')
        else:
            # FIXME: this is probably not very robust with regards to slightly non-matching meshes
            f0.set_allow_extrapolation(True)
            dolfin.project(f0, f.function_space(), function=f)


def setup_hooks(simulation):
    """
    Install the hooks that are given on the input file
    """
    simulation.log.info('Registering user-defined hooks')
    
    hooks = simulation.input.get_value('hooks', {}, 'dict(string:list)')
    
    def make_hook_from_code_string(name, code_string, description):
        runnable = RunnablePythonString(simulation, code_string, description)
        hook_data = simulation.io.get_persisted_dict('hook %r' % name)
        def hook(*args, **kwargs):
            runnable.run(hook_data=hook_data, **kwargs)
        return hook
    
    hook_types = [('pre_simulation', simulation.hooks.add_pre_simulation_hook),
                  ('post_simulation', simulation.hooks.add_post_simulation_hook),
                  ('pre_timestep', simulation.hooks.add_pre_timestep_hook),
                  ('post_timestep', simulation.hooks.add_post_timestep_hook),
                  ('matrix_ready', simulation.hooks.add_matrix_ready_hook)]
    names = set()
    
    for hook_name, register_hook in hook_types:
        for hook_info in hooks.get(hook_name, []):
            name = hook_info.get('name', 'unnamed')
            
            # Ensure that the hook name is unique among all hooks
            i, name2 = 0, name
            while name2 in names:
                i += 1
                name2 = name + str(i)
            name =  name2
            names.add(name)
                
            enabled = hook_info.get('enabled', True)
            description = '%s hook "%s"' % (hook_name, name)
            
            if not enabled:
                simulation.log.info('    Skipping disabled %s' % description)
                continue
            
            simulation.log.info('    Adding %s' % description)
            code_string = hook_info['code']
            hook = make_hook_from_code_string(name, code_string, description)
            register_hook(hook, 'User defined hook "%s"' % name)
            simulation.log.info('        ' + description)


def summarise_simulation_after_running(simulation, success):
    """
    Print a summary of the time spent on each part of the simulation
    """
    simulation.log.debug('\nGlobal simulation data at end of simulation:')
    for key, value in sorted(simulation.data.items()):
        simulation.log.debug('%20s = %s' % (key, repr(type(value))[:57]))
    
    # Add details on the time spent in each part of the simulation to the log
    clear = simulation.input.get_value('clear_timings_at_end', True, 'bool')
    log_timings(simulation, clear)
    
    # Show the total duration
    tottime = time.time() - simulation.t_start
    h = int(tottime/60**2)
    m = int((tottime - h*60**2)/60)
    s = tottime - h*60**2 - m*60
    humantime = '%d hours %d minutes and %d seconds' % (h, m, s)
    simulation.log.info('\nSimulation done in %.3f seconds (%s)' % (tottime, humantime))
    
    simulation.log.info("\nCurrent time: %s" % time.strftime('%Y-%m-%d %H:%M:%S'))
