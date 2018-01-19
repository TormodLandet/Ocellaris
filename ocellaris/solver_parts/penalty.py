import dolfin
from dolfin import cells, Constant, avg


def define_penalty(mesh, P, k_min, k_max, boost_factor=3, exponent=1):
    """
    Define the penalty parameter used in the Poisson equations
    
    Arguments:
        mesh: the mesh used in the simulation
        P: the polynomial degree of the unknown
        k_min: the minimum diffusion coefficient
        k_max: the maximum diffusion coefficient
        boost_factor: the penalty is multiplied by this factor
        exponent: set this to greater than 1 for superpenalisation
    """
    assert k_max >= k_min
    ndim = mesh.geometry().dim()
    
    # Calculate geometrical factor used in the penalty
    geom_fac = 0
    for cell in cells(mesh):
        vol = cell.volume()
        area = sum(cell.facet_area(i) for i in range(ndim + 1))
        gf = area/vol
        geom_fac = max(geom_fac, gf)
    geom_fac = dolfin.MPI.max(dolfin.MPI.comm_world, float(geom_fac))
    
    penalty = boost_factor * k_max**2/k_min * (P + 1)*(P + ndim)/ndim * geom_fac**exponent
    return penalty


def navier_stokes_stabilization_penalties(simulation, nu, velocity_continuity_factor_D12=0,
                                          pressure_continuity_factor=0, no_coeff=False):
    """
    Calculate the stabilization parameters needed in the DG scheme
    """
    ndim = simulation.ndim
    mpm = simulation.multi_phase_model
    mesh = simulation.data['mesh']
    
    if no_coeff:
        mu_min = mu_max = 1.0
    else:
        mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()
    
    P = simulation.data['Vu'].ufl_element().degree()
    penalty_dS = define_penalty(mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0)
    penalty_ds = penalty_dS * 2
    simulation.log.info('    DG SIP penalty:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds))
    
    if velocity_continuity_factor_D12:
        D12 = Constant([velocity_continuity_factor_D12]*ndim)
    else:
        D12 = Constant([0]*ndim)
    
    if pressure_continuity_factor:
        h = simulation.data['h']
        h = Constant(1.0)
        D11 = avg(h/nu)*Constant(pressure_continuity_factor)
    else:
        D11 = None
    
    return Constant(penalty_dS), Constant(penalty_ds), D11, D12
