# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import os
import sys
import subprocess
import pytest

MY_DIR = os.path.abspath(os.path.dirname(__file__))
SRC_DIR = os.path.join(MY_DIR, '..')


def pytest_generate_tests(metafunc):
    """
    Setup input files to test (all *.inp files in the same
    directory as this file and recursively in subdirctories)
    """

    def find_inp_files(dirname):
        path = os.path.abspath(os.path.join(SRC_DIR, dirname))
        for fn in os.listdir(path):
            pth_rel = os.path.join(dirname, fn)
            pth_abs = os.path.join(path, fn)
            if os.path.isdir(pth_abs):
                yield from find_inp_files(pth_rel)
            if fn.endswith('inp'):
                yield pth_rel[2:]

    inpfiles = list(find_inp_files('./demos'))
    metafunc.parametrize("demo_inp_file", inpfiles)


def test_setup_demo(demo_inp_file, monkeypatch):
    "Run setup on demo input file"

    # Skip some tests
    skip = ['linear_sloshing.inp', 'falling_sphere.inp']
    need_gmsh = ['cylinder.inp']
    slow = [
        'internal_soliton.inp',
        'dead_water_2D.inp',
        'dead_water_3D.inp',
        'wave_tank_vof_3D.inp',
    ]
    skip_patterns = skip + need_gmsh
    if os.environ.get('OCELLARIS_RUN_SLOW_TEST') != '1':
        skip_patterns += slow
    for pattern in skip_patterns:
        if pattern in demo_inp_file:
            return pytest.skip()

    # Run from the same directory as the demo file is located
    file_pth = os.path.abspath(demo_inp_file)
    dir_name = os.path.dirname(file_pth)
    base_name = os.path.basename(file_pth)
    monkeypatch.chdir(dir_name)

    # Some of the quicker simulations can run a couple of time steps
    tmax_for_test_run = {
        'lid_driven_cavity_flow.inp': 0.01,
        'taylor-green.inp': 0.1,
        'flow_around_ocellaris.inp': 0.001,
        'dam_break_2D.inp': 0.001,
        'dam_break_3D.inp': 0.0002,
    }
    tmax = tmax_for_test_run.get(base_name, 0.0)

    # Run as command (like the user would from the command line)
    # $ python3 -m ocellaris DEMOFILE --set-inp time/tmax=0
    interpreter = sys.executable
    mainfile = os.path.join(SRC_DIR, 'ocellaris', '__main__.py')
    cmd = [interpreter, mainfile, base_name]
    cmd.extend(['--set-inp', 'time/tmax=%r' % tmax])
    cmd.extend(['--set-inp', 'solver/num_inner_iter=2'])
    subprocess.check_call(cmd)

    # Run as script (can lead to inter-demo disturbances with dolfin globals etc)
    # from ocellaris import Simulation, setup_simulation
    # sim = Simulation()
    # sim.input.read_yaml(demo_inp_file)
    # success = setup_simulation(sim)
    # assert success
