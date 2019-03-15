# Copyright (C) 2017-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import dolfin
import pytest


@pytest.mark.parametrize('D', (2, 3))
def test_cell_midpoints(D):
    from ocellaris.solver_parts.slope_limiter.limiter_cpp_utils import SlopeLimiterInput

    if D == 2:
        mesh = dolfin.UnitSquareMesh(4, 4)
    else:
        mesh = dolfin.UnitCubeMesh(2, 2, 2)

    Vx = dolfin.FunctionSpace(mesh, 'DG', 2)
    V0 = dolfin.FunctionSpace(mesh, 'DG', 0)

    py_inp = SlopeLimiterInput(mesh, Vx, V0)
    cpp_inp = py_inp.cpp_obj

    all_ok = True
    for cell in dolfin.cells(mesh):
        cid = cell.index()
        mp = cell.midpoint()
        cpp_mp = cpp_inp.cell_midpoints[cid]
        for d in range(D):
            ok = dolfin.near(mp[d], cpp_mp[d])
            if not ok:
                print(
                    '%3d %d - %10.3e %10.3e' % (cid, d, mp[d], cpp_mp[d]),
                    '<--- ERROR' if not ok else '',
                )
                all_ok = False

    assert all_ok


def test_global_bounds():
    from ocellaris.solver_parts.slope_limiter.limiter_cpp_utils import SlopeLimiterInput

    mesh = dolfin.UnitSquareMesh(4, 4)

    Vx = dolfin.FunctionSpace(mesh, 'DG', 2)
    V0 = dolfin.FunctionSpace(mesh, 'DG', 0)

    py_inp = SlopeLimiterInput(mesh, Vx, V0)
    cpp_inp = py_inp.cpp_obj

    assert cpp_inp.global_min < -1e100
    assert cpp_inp.global_max > +1e100

    py_inp.set_global_bounds(0, 1)

    assert cpp_inp.global_min == 0
    assert cpp_inp.global_max == 1


def test_limit_cell_flag():
    from ocellaris.solver_parts.slope_limiter.limiter_cpp_utils import SlopeLimiterInput

    mesh = dolfin.UnitSquareMesh(4, 4)
    Ncells = mesh.topology().ghost_offset(2)

    Vx = dolfin.FunctionSpace(mesh, 'DG', 2)
    V0 = dolfin.FunctionSpace(mesh, 'DG', 0)

    py_inp = SlopeLimiterInput(mesh, Vx, V0)
    cpp_inp = py_inp.cpp_obj

    limit_cell = numpy.zeros(Ncells, dtype=numpy.intc)
    limit_cell[3] = 1
    py_inp.set_limit_cell(limit_cell)

    assert cpp_inp.limit_cell.shape == (Ncells,)
    assert cpp_inp.limit_cell[0] == 0
    assert cpp_inp.limit_cell[3] == 1
