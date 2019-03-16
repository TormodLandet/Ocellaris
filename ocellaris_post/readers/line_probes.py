# Copyright (C) 2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import numpy


def read_line_probes(res):
    inp = res.input
    res.line_probes = {}
    if 'probes' not in inp:
        return

    for probe in inp['probes']:
        if not (probe.get('enabled', True) and probe.get('type', '') == 'LineProbe'):
            continue
        name = probe['name']
        field_name = probe['field']
        file_name_postfix = '_' + probe['file_name']
        file_name = res.get_file_path(file_name_postfix)
        target_name = probe.get('target_name', 'Target')
        target_abcissa = probe.get('target_abcissa', [])
        target_ordinate = probe.get('target_ordinate', [])
        probe = LineProbe(name, field_name, file_name, target_name, target_abcissa, target_ordinate)
        res.line_probes[name] = probe


class LineProbe(object):
    def __init__(self, name, field_name, file_name, target_name, target_abcissa, target_ordinate):
        self.name = name
        self.field_name = field_name
        self.file_name = file_name
        self._cache = None

        self.target_name = target_name
        self.target_abcissa = target_abcissa
        self.target_ordinate = target_ordinate

    def reload(self):
        self._cache = None

    def get_line_probes(self, cache=True):
        if cache and self._cache is not None:
            return self._cache

        timesteps = []
        data = []

        with open(self.file_name, 'rt') as f:
            description = f.readline()[1:].strip()
            xvec = numpy.array([float(v) for v in f.readline()[5:].split()])
            yvec = numpy.array([float(v) for v in f.readline()[5:].split()])
            zvec = numpy.array([float(v) for v in f.readline()[5:].split()])
            assert xvec.shape == yvec.shape == zvec.shape

            # Separator line
            line = f.readline()

            while line:
                line = f.readline()
                wds = line.split()
                if len(wds) != len(xvec) + 1:
                    # Missing data on this line, file not fully written?
                    continue

                timesteps.append(float(wds[0]))
                values = numpy.array([float(v) for v in wds[1:]])
                data.append(values)

        res = (description, numpy.array(timesteps), xvec, yvec, zvec, data)
        if cache:
            self._cache = res
        return res
