import time
import json
import dolfin
from .multiphase import get_multi_phase_model
from .solvers import get_solver
from .boundary_conditions import BoundaryRegion

def run_simulation(simulation, verbose=True):
    """
    Prepare and run a simulation
    """
    if verbose:
        print 'Preparing simulation ...'
        print
        t1 = time.time()
    
    # Load the mesh
    load_mesh(simulation)
    
    # Create function spaces
    make_function_spaces(simulation)
    
    # Load the boundary conditions
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
    if verbose:
        print 'Preparing simulation done in %.3f seconds' % (time.time() - t1)
        print 
        print 'Simulation configuration:'
        print '{:-^40}'.format(' input begin ')
        print json.dumps(simulation.input, indent=4)
        print '{:-^40}'.format(' input end ')
        print
        print "Running simulation ..."
        t1 = time.time()

    # Run the simulation 
    solver.run()
    
    if simulation.input.get('debug', True):
        print '\nGlobal simulation data at end of simulation:'
        for key, value in sorted(simulation.data.items()):
            print '%20s = %s' % (key, repr(type(value))[:57])
        print
    
    if verbose:
        print 'Simulation done in %.3f seconds' % (time.time() - t1)

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
    for index, area in enumerate(simulation.input['boundary_conditions']):
        part = BoundaryRegion(simulation, marker, index)
        boundary.append(part)
    simulation.data['boundary'] = boundary
    
    # Create a boundary meassure that is aware of the marked regions
    ds = dolfin.Measure('ds')[marker]
    simulation.data['ds'] = ds
    
    print 'Dirichlet boundary conditions:'
    for key, vals in sorted(simulation.data['dirichlet_bcs'].items()):
        print '  ', key
        for v in vals:
            print '    ', v 
    exit()