# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy


def read_surfaces(res):
    inp = res.input
    res.surfaces = {}
    if 'probes' not in inp:
        return

    for probe in inp['probes']:
        if not (probe.get('enabled', True) and probe.get('type', '') == 'IsoSurface'):
            continue
        name = probe['name']
        field_name = probe['field']
        value = probe['value']
        file_name_postfix = probe['file_name']
        file_name = res.get_file_path(file_name_postfix)
        isosurf = IsoSurfaces(name, field_name, value, file_name)
        res.surfaces[name] = isosurf


class IsoSurfaces(object):
    def __init__(self, name, field_name, value, file_name):
        self.name = name
        self.field_name = field_name
        self.value = value
        self.file_name = file_name
        self._cache = None

    def reload(self):
        self._cache = None

    def get_surfaces(self, cache=True):
        if cache and self._cache is not None:
            return self._cache

        timesteps = []
        data = []

        with open(self.file_name, 'rt') as f:
            description = f.readline()[1:].strip()
            value = float(f.readline().split()[-1])
            dim = int(f.readline().split()[-1])

            line = f.readline()
            while line:
                wds = line.split()
                try:
                    time = float(wds[1])
                    nsurf = int(wds[3])
                except Exception:
                    break

                if nsurf == 0:
                    timesteps.append(time)
                    data.append([])
                    line = f.readline()
                    continue

                datalines = [f.readline() for _ in range(nsurf * 3)]
                if not datalines[-1]:
                    break
                timesteps.append(time)
                data.append([])
                for i in range(nsurf):
                    xvals = [float(v) for v in datalines[i * 3 + 0].split()]
                    yvals = [float(v) for v in datalines[i * 3 + 1].split()]
                    zvals = [float(v) for v in datalines[i * 3 + 2].split()]
                    data[-1].append((xvals, yvals, zvals))

                line = f.readline()

        res = (description, value, dim, numpy.array(timesteps), data)
        if cache:
            self._cache = res
        return res
