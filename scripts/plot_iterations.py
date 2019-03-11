# Copyright (C) 2016-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0
"""
Plot iteration convergence logs from an Ocellaris restart file 
"""
import h5py
from matplotlib import pyplot


def read_iterations(h5_file_name, derived=True):
    hdf = h5py.File(h5_file_name, 'r')

    meta = hdf['/ocellaris']

    log = []
    i = 0
    while True:
        logname = 'full_log_%d' % i
        i += 1
        if not logname in meta.attrs:
            break
        log.append(meta.attrs[logname])

    log = ''.join(log).split('\n')
    errors_u = []
    errors_p = []
    errors_div = []
    iterations_per_timestep = []
    it = 0
    for line in log:
        if 'Reports for timestep' in line:
            iterations_per_timestep.append(it)
            it = 0
        elif 'iteration' in line and 'err u*' in line:
            it += 1
            wds = line.split()
            errors_u.append(float(wds[6]))
            errors_p.append(float(wds[10]))
            errors_div.append(float(wds[-1]))

    if it != 0:
        iterations_per_timestep.append(it)

    return errors_u, errors_p, errors_div, iterations_per_timestep


def plot_iterations(file_name):
    errors_u, errors_p, errors_div, iterations_per_timestep = read_iterations(file_name)

    fig = pyplot.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax1.plot(errors_u)
    ax2.plot(errors_p)
    ax3.plot(errors_div)
    for split in iterations_per_timestep:
        ax1.axvline(split, color='r')
        ax2.axvline(split, color='r')
        ax3.axvline(split, color='r')

    # ax1.set_ylim(0, 0.1)
    # ax2.set_ylim(0, 0.1)
    # ax3.set_ylim(0, 0.1)

    fig.tight_layout()
    pyplot.show()


if __name__ == '__main__':
    import sys

    h5_file_name = sys.argv[1]
    plot_iterations(h5_file_name)
    pyplot.show()
