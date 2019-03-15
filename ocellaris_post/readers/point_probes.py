# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy
import os


def read_point_probes(res):
    inp = res.input
    res.point_probes = {}
    if 'probes' not in inp:
        return

    for probe in inp['probes']:
        if not (probe.get('enabled', True) and probe.get('type', '') == 'PointProbe'):
            continue
        name = probe['name']
        prefix = res.get_file_path('', check=False)
        file_name = probe.get('file_name', prefix + '_%s.tsv' % name)

        if not os.path.isfile(file_name):
            res.warnings.append('PointProbe file not found: %s' % file_name)
            continue

        probes = PointProbes(name, file_name)
        res.point_probes[name] = probes


class PointProbes(object):
    def __init__(self, name, file_name):
        self.name = name
        self.file_name = file_name
        self.reload()

    def reload(self):
        self._cache = None
        with open(self.file_name, 'rt') as tsv:
            header = tsv.readline()
        self.probe_names = header.strip().split('\t')[1:]

    def get_probe(self, probe_name, cache=True):
        if cache and self._cache is not None:
            timesteps, probes = self._cache
            return timesteps, probes[probe_name]

        data = numpy.loadtxt(self.file_name, skiprows=1).T
        timesteps = data[0]
        probes = {}
        for i, pname in enumerate(self.probe_names):
            probes[pname] = data[i + 1]

        if cache:
            self._cache = timesteps, probes
        return timesteps, probes[probe_name]
