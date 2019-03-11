# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from dolfin import cells, Constant, avg


def define_penalty(mesh, P, k_min, k_max, boost_factor=3, exponent=1):
    """
    Define the penalty parameter used in the Poisson equations

    Spatially constant version, returns a constant expression

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
        gf = area / vol
        geom_fac = max(geom_fac, gf)
    geom_fac = dolfin.MPI.max(dolfin.MPI.comm_world, float(geom_fac))

    penalty = (
        boost_factor
        * k_max ** 2
        / k_min
        * (P + 1)
        * (P + ndim)
        / ndim
        * geom_fac ** exponent
    )
    return penalty


def define_spatially_varying_penalty(
    simulation, P, k_min, k_max, boost_factor=3, exponent=1
):
    """
    Define the penalty parameter used in the Poisson equations

    Spatially varying version, returns a DGT0 function

    Arguments:
        mesh: the mesh used in the simulation
        P: the polynomial degree of the unknown
        k_min: the minimum diffusion coefficient
        k_max: the maximum diffusion coefficient
        boost_factor: the penalty is multiplied by this factor
        exponent: set this to greater than 1 for superpenalisation
    """
    assert k_max >= k_min
    mesh = simulation.data['mesh']
    ndim = mesh.geometry().dim()

    # Compute the constant part of the penalty
    pconst = boost_factor * k_max ** 2 / k_min * (P + 1) * (P + ndim) / ndim

    # Compute the spatially varying penalty
    V = dolfin.FunctionSpace(mesh, 'DG', 0)
    dm = V.dofmap()
    penalty_func = dolfin.Function(V)
    arr = penalty_func.vector().get_local()
    for cell in dolfin.cells(mesh):
        vol = cell.volume()
        area = sum(cell.facet_area(i) for i in range(ndim + 1))
        geom_fac = area / vol
        dof, = dm.cell_dofs(cell.index())
        arr[dof] = pconst * geom_fac ** exponent

    penalty_func.vector().set_local(arr)
    penalty_func.vector().apply('insert')

    # Optionally plot the penalty function to file
    if simulation.input.get_value('output/plot_elliptic_penalty', False, 'bool'):
        prefix = simulation.input.get_value('output/prefix', '', 'string')
        pfile = prefix + '_elliptic_dg_penalty.xdmf'
        simulation.log.info('    Plotting elliptic DG penalty to XDMF file %r' % pfile)
        penalty_func.rename('penalty', 'penalty')
        with dolfin.XDMFFile(mesh.mpi_comm(), pfile) as xdmf:
            xdmf.write(penalty_func)

    return penalty_func


def navier_stokes_stabilization_penalties(
    simulation,
    nu,
    velocity_continuity_factor_D12=0,
    pressure_continuity_factor=0,
    no_coeff=False,
):
    """
    Calculate the stabilization parameters needed in the DG scheme
    """
    ndim = simulation.ndim
    mpm = simulation.multi_phase_model
    mesh = simulation.data['mesh']
    use_const = simulation.input.get_value(
        'solver/spatially_constant_penalty', False, 'bool'
    )

    if no_coeff:
        mu_min = mu_max = 1.0
    else:
        mu_min, mu_max = mpm.get_laminar_dynamic_viscosity_range()

    P = simulation.data['Vu'].ufl_element().degree()
    if use_const:
        simulation.log.info('    Using spatially constant elliptic penalty')
        penalty_dS = define_penalty(
            mesh, P, mu_min, mu_max, boost_factor=3, exponent=1.0
        )
        penalty_ds = penalty_dS * 2
        simulation.log.info(
            '    DG SIP penalty:  dS %.1f  ds %.1f' % (penalty_dS, penalty_ds)
        )
        penalty_dS = Constant(penalty_dS)
        penalty_ds = Constant(penalty_ds)

    else:
        simulation.log.info('    Using spatially varying elliptic penalties')
        penalty_dS = define_spatially_varying_penalty(
            simulation, P, mu_min, mu_max, boost_factor=3, exponent=1.0
        )
        penalty_ds = penalty_dS * 2
        penalty_dS = dolfin.conditional(
            dolfin.lt(penalty_dS('+'), penalty_dS('-')),
            penalty_dS('-'),
            penalty_dS('+'),
        )

    if velocity_continuity_factor_D12:
        D12 = Constant([velocity_continuity_factor_D12] * ndim)
    else:
        D12 = Constant([0] * ndim)

    if pressure_continuity_factor:
        h = simulation.data['h']
        h = Constant(1.0)
        D11 = avg(h / nu) * Constant(pressure_continuity_factor)
    else:
        D11 = None

    return penalty_dS, penalty_ds, D11, D12
