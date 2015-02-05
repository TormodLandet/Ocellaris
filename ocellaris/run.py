import dolfin
from .multiphase import get_multi_phase_model
from .solvers import get_solver

def run_simulation(simulation):
    """
    Prepare a simulation for running and then run it
    """
    # Load the mesh
    load_mesh(simulation)
    
    # Get the density and viscosity properties from the multi phase model
    multiphase_model_name = simulation.input.get('multiphase_model', 'SinglePhase')
    multiphase_model = get_multi_phase_model(multiphase_model_name)(simulation)
    simulation.data['rho'] = multiphase_model.get_density()
    simulation.data['nu'] = multiphase_model.get_laminar_kinematic_viscosity()
    
    # Get the solver
    solver_name = simulation.input['solver']['type']
    solver = get_solver(solver_name)(simulation)

    # Run the simulation
    solver.run()

def load_mesh(simulation):
    """
    Get the mesh from the input file
    """
    mesh_type = simulation.input['mesh']['type']
    
    if mesh_type == 'Rectangle':
        startx = simulation.input['mesh'].get('startx', 0)
        starty = simulation.input['mesh'].get('starty', 0)
        endx = simulation.input['mesh'].get('endx', 1)
        endy = simulation.input['mesh'].get('endy', 1)
        Nx = int(simulation.input['mesh']['Nx'])
        Ny = int(simulation.input['mesh']['Ny'])
        diagonal = simulation.input['mesh'].get('diagonal', 'left/right')
        mesh = dolfin.RectangleMesh(startx, starty, endx, endy, Nx, Ny, diagonal)
    
    simulation.set_mesh(mesh)
