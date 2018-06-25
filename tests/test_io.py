"""
Test restart file handling and general IO
"""
import numpy
import os
from ocellaris import Simulation, setup_simulation
from poisson_solver import BASE_INPUT as BASE_INPUT_PHI
from helpers import mpi_tmpdir
import pytest


BASE_INPUT_VELPRES = """
ocellaris:
    type: input
    version: 1.0

user_code:
    constants:
        N: 3
        a: 3.1415926535897932384626433
        d: 3.1415926535897932384626433
        u0a: '-a * (exp(a * x[0]) * sin(a * x[1] + d * x[2]) + exp(a * x[2]) * cos(a * x[0] + d * x[1])) * exp(-d * d * t)'
        u1a: '-a * (exp(a * x[1]) * sin(a * x[2] + d * x[0]) + exp(a * x[0]) * cos(a * x[1] + d * x[2])) * exp(-d * d * t)'
        u2a: '-a * (exp(a * x[2]) * sin(a * x[0] + d * x[1]) + exp(a * x[1]) * cos(a * x[2] + d * x[0])) * exp(-d * d * t)'
        pa:  |
          -a * a / 2 * (
              exp(2 * a * x[0]) + exp(2 * a * x[1]) + exp(2 * a * x[2]) +
              2 * sin(a * x[0] + d * x[1]) * cos(a * x[2] + d * x[0]) * exp(a * (x[1] + x[2])) + 
              2 * sin(a * x[1] + d * x[2]) * cos(a * x[0] + d * x[1]) * exp(a * (x[2] + x[0])) + 
              2 * sin(a * x[2] + d * x[0]) * cos(a * x[1] + d * x[2]) * exp(a * (x[0] + x[1]))
          ) * exp(-2 * d * d * t)

mesh:
    type: Box
    Nx: py$ N
    Ny: py$ N
    Nz: py$ N

solver:
    type: AnalyticalSolution

output:
    log_enabled: yes
    solution_properties: off
    save_restart_file_at_end: off

initial_conditions:
    up0:
        cpp_code: py$ u0a
    up1:
        cpp_code: py$ u1a
    up2:
        cpp_code: py$ u2a
    p:
        cpp_code: py$ pa

# Dummy values below here are just to make Ocellaris run. The
# solution framework has some assumptions that all solvers
# are (quasi) time stepping flow solvers
time:
    dt: 1.0
physical_properties:
    nu: 1.0
boundary_conditions: []
"""


def test_restart_file_io(tmpdir_factory):
    dir_name = mpi_tmpdir(tmpdir_factory, 'test_restart_file_io')
    prefix = os.path.join(dir_name, 'ocellaris')

    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT_PHI)
    sim.input.set_value('output/prefix', prefix)
    sim.input.set_value('time/tstart', 42.0)
    setup_simulation(sim)

    # Fill in the phi function
    phi = sim.data['phi']
    phi_arr = phi.vector().get_local()
    phi_arr[:] = numpy.random.rand(*phi_arr.shape)
    phi.vector().set_local(phi_arr)
    phi.vector().apply('insert')

    # Save restart file
    file_name = sim.io.write_restart_file()
    assert file_name.startswith(prefix)

    # Load input from restart file
    sim2 = Simulation()
    sim2.io.load_restart_file_input(file_name)
    assert sim2.input.get_value('time/tstart') == 42.0
    assert str(sim.input) == str(sim2.input)

    # Load phi from restart file
    setup_simulation(sim2)
    sim2.io.load_restart_file_results(file_name)
    phi2 = sim2.data['phi']
    phi2_arr = phi2.vector().get_local()

    # FIXME: make tests work in parallel
    if sim2.data['mesh'].mpi_comm().size == 1:
        assert all(phi_arr == phi2_arr)
        assert sim.data['mesh'].hash() == sim2.data['mesh'].hash()


@pytest.mark.parametrize("iotype", ['vtk', 'xdmf'])
def test_plot_io_3D(iotype, tmpdir_factory):
    dir_name = mpi_tmpdir(tmpdir_factory, 'test_plot_io_3D')
    prefix = os.path.join(dir_name, 'ocellaris')
    N = 4

    sim = Simulation()
    sim.input.read_yaml(yaml_string=BASE_INPUT_VELPRES)
    sim.input.set_value('output/prefix', prefix)
    sim.input.set_value('user_code/constants/N', N)
    setup_simulation(sim)
    sim.io.setup()  # Normally called in sim._at_start_of_simulation()

    for d in range(3):
        sim.data['u%d' % d].assign(sim.data['up%d' % d])

    ncell = N ** 3 * 6
    npoint = ncell * 10

    if iotype == 'vtk':
        file_name = sim.io.lvtk.write()
        assert file_name.startswith(prefix)
        assert file_name.endswith('.vtk')
        assert os.path.isfile(file_name)

        with open(file_name, 'rt') as f:
            for line in f:
                if line.startswith('POINTS'):
                    pline = line
                elif line.startswith('CELLS'):
                    cline = line
                elif line.startswith('CELL_TYPES'):
                    tline = line
                elif line.startswith('POINT_DATA'):
                    dline = line
        assert pline.strip() == 'POINTS %d float' % npoint
        assert cline.strip() == 'CELLS %d %d' % (ncell, ncell * 11)
        assert tline.strip() == 'CELL_TYPES %d' % ncell
        assert dline.strip() == 'POINT_DATA %d' % (ncell * 10)

    elif iotype == 'xdmf':
        file_name = sim.io.xdmf.write()
        assert file_name.startswith(prefix)
        assert file_name.endswith('.xdmf')
        assert os.path.isfile(file_name)
