# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import sys
import os
import time
import importlib
import dolfin
from ocellaris.solvers import get_solver
from ocellaris.probes import setup_probes
from ocellaris.utils import (
    interactive_console_hook,
    ocellaris_error,
    ocellaris_interpolate,
    log_timings,
    verify_field_variable_definition,
    load_meshio_mesh,
    build_distributed_mesh,
)
from ocellaris.utils import RunnablePythonString, OcellarisCppExpression
from ocellaris.utils import verify_key
from ocellaris.solver_parts import (
    BoundaryRegion,
    get_multi_phase_model,
    get_known_field,
    MeshMorpher,
    add_forcing_zone,
)


def setup_simulation(simulation):
    """
    Prepare an Ocellaris simulation for running
    """
    # Force all log output to be flushed straight away (for better debugging)
    force_flush_old = simulation.log.force_flush_all
    simulation.log.force_flush_all = True

    # Setup user code (custom sys.path, custom module imports)
    setup_user_code(simulation)

    simulation.log.info('Preparing simulation')
    simulation.log.info(
        'Output prefix is: %s' % simulation.input.get_value('output/prefix', '', 'string')
    )
    simulation.log.info("Current time: %s\n" % time.strftime('%Y-%m-%d %H:%M:%S'))
    t_start = time.time()

    # Set linear algebra backend to PETSc
    setup_fenics(simulation)

    # Make time and timestep available in expressions for the initial conditions etc
    simulation.log.info('Creating time simulation')
    simulation.time = simulation.input.get_value('time/tstart', 0.0, 'float')
    simulation.dt = simulation.input.get_value('time/dt', None, 'float')
    simulation.dt_prev = simulation.dt
    if simulation.dt is None:
        simulation.dt = 1e100
        simulation.log.warning(
            'No timestep (time/dt) given. Using %r' % simulation.dt
            + ', make sure to change this if this is not a static simulation!'
        )
    assert simulation.dt > 0

    # Preliminaries, before setup begins

    # Get the multi phase model class
    multiphase_model_name = simulation.input.get_value(
        'multiphase_solver/type', 'SinglePhase', 'string'
    )
    multiphase_class = get_multi_phase_model(multiphase_model_name)

    # Get the main solver class
    solver_name = simulation.input.get_value('solver/type', required_type='string')
    simulation.log.info('Using equation solver %s' % solver_name)
    solver_class = get_solver(solver_name)

    ###########################################################################
    # Setup the Ocellaris simulation

    if not simulation.restarted:
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
    setup_function_spaces(simulation, solver_class, multiphase_class)

    # Load the mesh morpher used for prescribed mesh velocities and ALE multiphase solvers
    simulation.mesh_morpher = MeshMorpher(simulation)

    # Setup physical constants and multi-phase model (g, rho, nu, mu)
    setup_physical_properties(simulation, multiphase_class)

    # Setup known fields (incomming waves etc)
    setup_known_fields(simulation)

    # Load the boundary conditions. This must be done
    # before creating the solver as the solver needs
    # the Neumann conditions to define weak forms
    setup_boundary_conditions(simulation)

    # Create momentum sources (useful for MMS tests etc)
    setup_sources(simulation)

    # Create the solver
    simulation.log.info('Initialising %s solver' % solver_name)
    simulation.solver = solver_class(simulation)

    # Setup postprocessing probes
    setup_probes(simulation)

    # Initialise the fields
    if not simulation.restarted:
        setup_initial_conditions(simulation)

    # Setup the solution properties
    simulation.solution_properties.setup()

    # Setup any hooks that may be present on the input file
    setup_hooks(simulation)

    # Setup the interactive console to optionally run at the end of each timestep
    simulation.hooks.add_post_timestep_hook(
        lambda: interactive_console_hook(simulation), 'Interactive console commands'
    )

    # Setup the summary to show after the simulation
    def hook(success):
        return summarise_simulation_after_running(simulation, success)

    simulation.hooks.add_post_simulation_hook(hook, 'Summarise simulation')

    # Show time spent setting up the solver
    simulation.log.info('\nPreparing simulation done in %.3f seconds' % (time.time() - t_start))

    # Show all registered hooks
    simulation.hooks.show_hook_info()
    simulation.log.force_flush_all = force_flush_old


def setup_user_code(simulation):
    """
    Setup custom path and user code imports
    """
    # Add the directory where the input file resides to the path
    if simulation.input.file_name is not None:
        fn = os.path.abspath(simulation.input.file_name)
        inp_dir = os.path.dirname(fn)
        sys.path.insert(0, inp_dir)

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
        ocellaris_error(
            'Missing PETSc',
            'DOLFIN has not been configured with PETSc ' 'which is needed by Ocellaris.',
        )
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
        simulation.log.warning('(deprecated, please use meshio reader)')

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

    elif mesh_type == 'XDMF':
        simulation.log.info('Creating mesh from XDMF file')
        simulation.log.warning('(deprecated, please use meshio reader)')

        mesh_file = inp.get_value('mesh/mesh_file', required_type='string')

        # Load the mesh from file
        pth = inp.get_input_file_path(mesh_file)
        mesh = dolfin.Mesh(comm)
        with dolfin.XDMFFile(comm, pth) as xdmf:
            xdmf.read(mesh, False)

    elif mesh_type == 'HDF5':
        simulation.log.info('Creating mesh from DOLFIN HDF5 file')
        h5_file_name = inp.get_value('mesh/mesh_file', required_type='string')

        with dolfin.HDF5File(comm, h5_file_name, 'r') as h5:
            # Read mesh
            mesh = dolfin.Mesh(comm)
            h5.read(mesh, '/mesh', False)

            # Read facet regions
            if h5.has_dataset('/mesh_facet_regions'):
                mesh_facet_regions = dolfin.FacetFunction('size_t', mesh)
                h5.read(mesh_facet_regions, '/mesh_facet_regions')
            else:
                mesh_facet_regions = None

    elif mesh_type == 'meshio':
        simulation.log.info('Creating mesh with meshio reader')
        file_name = inp.get_value('mesh/mesh_file', required_type='string')
        file_type = inp.get_value('mesh/meshio_type', None, required_type='string')
        sort_order = inp.get_value('mesh/sort_order', None, required_type='list(int)')

        if sort_order:
            simulation.log.info(
                '    Ordering mesh elements by ' + ', then '.join('xyz'[i] for i in sort_order)
            )

        # Read mesh on rank 0
        t1 = time.time()
        mesh = dolfin.Mesh(comm)
        if comm.rank == 0:
            # Read a mesh file by use of meshio
            physical_regions = load_meshio_mesh(mesh, file_name, file_type, sort_order)
            simulation.log.info(
                '    Read mesh with %d cells in %.2f seconds' % (mesh.num_cells(), time.time() - t1)
            )
        else:
            physical_regions = None

        # Distribute the mesh
        if comm.size > 1:
            physical_regions = comm.bcast(physical_regions)
            t1 = time.time()
            simulation.log.info('    Distributing mesh to %d processors' % comm.size)
            build_distributed_mesh(mesh)
            simulation.log.info('    Distributed mesh in %.2f seconds' % (time.time() - t1))
    else:
        ocellaris_error('Unknown mesh type', 'Mesh type %r is not supported' % mesh_type)

    # Optionally move the mesh (for simple grading etc)
    move = inp.get_value('mesh/move', None, required_type='list(string)')
    if move is not None:
        simulation.log.info('    Moving mesh')
        if len(move) != mesh.geometry().dim():
            ocellaris_error(
                'Mesh move not correct',
                'Length of move field is %d while geometric dimension is %r'
                % (len(move), mesh.geometry().dim()),
            )
        e_move = dolfin.Expression(move, degree=1)
        V_move = dolfin.VectorFunctionSpace(mesh, 'CG', 1)
        f_move = dolfin.interpolate(e_move, V_move)
        dolfin.ALE.move(mesh, f_move)
        mesh.bounding_box_tree().build(mesh)

    # Update the simulation
    simulation.set_mesh(mesh, mesh_facet_regions)

    # Load meshio facet regions
    if mesh_type == 'meshio' and physical_regions:
        mfr = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
        conn_FV = simulation.data['connectivity_FV']
        vcoords = mesh.coordinates()
        for ifacet in range(mfr.size()):
            facet_vertices = conn_FV(ifacet)
            key = sorted([tuple(vcoords[vidx]) for vidx in facet_vertices])
            number = physical_regions.get(tuple(key), 0)
            mfr[ifacet] = number
        # Store the loaded regions
        simulation.data['mesh_facet_regions'] = mfr
        mesh_facet_regions = mfr

    # Optionally plot mesh right after loading (for debugging)
    if simulation.input.get_value('output/plot_mesh', False, 'bool'):
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        pfile = prefix + '_mesh.xdmf'
        simulation.log.info('    Plotting mesh with current MPI ranks to XDMF file %r' % pfile)
        V0 = dolfin.FunctionSpace(mesh, 'DG', 0)
        ranks = dolfin.Function(V0)
        ranks.vector().set_local(ranks.vector().get_local() * 0 + comm.rank)
        ranks.vector().apply('insert')
        ranks.rename('MPI_rank', 'MPI_rank')
        with dolfin.XDMFFile(comm, pfile) as xdmf:
            xdmf.write(ranks)

    # Optionally plot facet regions to file
    if simulation.input.get_value('output/plot_facet_regions', False, 'bool'):
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        pfile = prefix + '_input_facet_regions.xdmf'
        simulation.log.info('    Plotting input mesh facet regions to XDMF file %r' % pfile)
        if mesh_facet_regions is None:
            simulation.log.warning('Cannot plot mesh facet regions, no regions found!')
        else:
            with dolfin.XDMFFile(comm, pfile) as xdmf:
                xdmf.write(mesh_facet_regions)


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
    for index, _ in enumerate(simulation.input.get_value('boundary_conditions', [], 'list(dict)')):
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
        length = dolfin.assemble(one * ds(region.mark_id, domain=mesh))
        pf = simulation.log.info if length > 0.0 else simulation.log.warning
        pf('    Boundary region %s has size %f' % (region.name, length))
    length0 = dolfin.assemble(one * ds(0, domain=mesh))
    pf = simulation.log.info if length0 == 0.0 else simulation.log.warning
    pf('    Boundary region UNMARKED has size %f' % length0)

    # Optionally plot boundary regions to file
    if simulation.input.get_value('output/plot_bcs', False, 'bool'):
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        pfile = prefix + '_boundary_regions.xdmf'
        simulation.log.info('    Plotting boundary regions to ' 'XDMF file %r' % pfile)
        with dolfin.XDMFFile(mesh.mpi_comm(), pfile) as xdmf:
            xdmf.write(marker)


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
            simulation.log.info(
                '    Function space %s has dimension %d (%s degree %d)'
                % (name, V.dim(), family, degree)
            )


def setup_physical_properties(simulation, multiphase_class):
    """
    Gravity vector and rho/nu/mu fields are created here
    """
    ndim = simulation.ndim
    g = simulation.input.get_value(
        'physical_properties/g',
        [0] * ndim,
        required_type='list(float)',
        required_length=simulation.ndim,
    )
    simulation.data['g'] = dolfin.Constant(g)

    # Get the density and viscosity properties from the multi phase model
    simulation.multi_phase_model = multiphase_class(simulation)

    simulation.data['rho'] = simulation.multi_phase_model.get_density(0)
    simulation.data['nu'] = simulation.multi_phase_model.get_laminar_kinematic_viscosity(0)
    simulation.data['mu'] = simulation.multi_phase_model.get_laminar_dynamic_viscosity(0)


def setup_known_fields(simulation):
    """
    Setup known fields like incoming waves etc
    """
    simulation.log.info('Creating known fields')
    Nfields = len(simulation.input.get_value('fields', [], 'list(dict)'))
    for i in range(Nfields):
        field_inp = simulation.input.get_value('fields/%d' % i, required_type='Input')
        name = field_inp.get_value('name', required_type='string')
        field_type = field_inp.get_value('type', required_type='string')
        if name in simulation.fields:
            ocellaris_error(
                'Field %s is defined multiple times' % name,
                'Input file contains multiple fields with the same name',
            )
        field_class = get_known_field(field_type)
        simulation.fields[name] = field_class(simulation, field_inp)


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
    Setup the momentum sources and forcing zones
    """
    # Standard momentum sources
    ms = simulation.input.get_value('momentum_sources', [], 'list(dict)')
    simulation.data['momentum_sources'] = sources = []
    if ms:
        simulation.log.info('Creating sources')
        for i in range(len(ms)):
            inp = simulation.input.get_value('momentum_sources/%d', required_type='Input')
            name = inp.get_value('name', required_type='string')
            degree = inp.get_value('degree', required_type='int')
            cpp_code = inp.get_value('cpp_code', required_type='list(string)')
            description = 'momentum source %r' % name
            simulation.log.info('    C++ %s' % description)

            expr = OcellarisCppExpression(simulation, cpp_code, description, degree, update=True)
            sources.append(expr)

    # Penalty forcing zones
    fz = simulation.input.get_value('forcing_zones', [], 'list(dict)')
    simulation.data['forcing_zones'] = fzones = {}
    if fz:
        simulation.log.info('Creating forcing zones')
        for i in range(len(fz)):
            inp = simulation.input.get_value('forcing_zones/%d' % i, required_type='Input')
            name = inp.get_value('name', required_type='string')
            ztype = inp.get_value('type', required_type='string')
            description = '%s zone %r' % (ztype, name)
            simulation.log.info('    %s' % description)
            add_forcing_zone(simulation, fzones, inp)


def setup_initial_conditions(simulation):
    """
    Setup the initial values for the fields

    NOTE: this is never run on simulation restarts!
    """
    simulation.log.info('Creating initial conditions')

    ic = simulation.input.get_value('initial_conditions', {}, 'Input')
    has_file = False
    for name in ic:
        name = str(name)

        if name == 'file':
            has_file = True
            continue
        elif 'p' not in name:
            ocellaris_error(
                'Invalid initial condition',
                'You have given initial conditions for %r but this does '
                'not seem to be a previous or pressure field.\n\n'
                'Valid names: up0, up1, ... , p, cp, rho_p, ...' % name,
            )
        elif name not in simulation.data:
            ocellaris_error(
                'Invalid initial condition',
                'You have given initial conditions for %r but this does '
                'not seem to be an existing field.' % name,
            )

        func = simulation.data[name]
        V = func.function_space()
        description = 'initial conditions for %r' % name

        if 'cpp_code' in ic[name]:
            cpp_code = ic.get_value('%s/cpp_code' % name, required_type='string!')
            simulation.log.info('    C++ %s' % description)
            ocellaris_interpolate(simulation, cpp_code, description, V, func)
        elif 'function' in ic[name]:
            vardef = ic.get_value('%s/function' % name, required_type='string!')
            simulation.log.info('    Field function %s' % description)
            f = verify_field_variable_definition(simulation, vardef, description)
            dolfin.project(f, V, function=func)
        else:
            ocellaris_error(
                'Invalid initial condition',
                'You have not given "cpp_code" or "function" for %r' % name,
            )

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
            simulation.log.warning(
                '    Found initial condition for %r in h5 file, but '
                'this function is not present in the simulation' % name
            )
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

    hook_types = [
        ('pre_simulation', simulation.hooks.add_pre_simulation_hook),
        ('post_simulation', simulation.hooks.add_post_simulation_hook),
        ('pre_timestep', simulation.hooks.add_pre_timestep_hook),
        ('post_timestep', simulation.hooks.add_post_timestep_hook),
        ('matrix_ready', simulation.hooks.add_matrix_ready_hook),
    ]
    names = set()

    for hook_name, register_hook in hook_types:
        for hook_info in hooks.get(hook_name, []):
            name = hook_info.get('name', 'unnamed')

            # Ensure that the hook name is unique among all hooks
            i, name2 = 0, name
            while name2 in names:
                i += 1
                name2 = name + str(i)
            name = name2
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
    h = int(tottime / 60 ** 2)
    m = int((tottime - h * 60 ** 2) / 60)
    s = tottime - h * 60 ** 2 - m * 60
    humantime = '%d hours %d minutes and %d seconds' % (h, m, s)
    simulation.log.info('\nSimulation done in %.3f seconds (%s)' % (tottime, humantime))

    simulation.log.info("\nCurrent time: %s" % time.strftime('%Y-%m-%d %H:%M:%S'))
