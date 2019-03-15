# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin


class DebugIO:
    def __init__(self, simulation):
        """
        Debugging IO routines. These routines are not used during normal runs
        (and are in general rarely used/tested at all)
        """
        self.simulation = simulation

    def save_la_objects(self, file_name, **kwargs):
        """
        Save a dictionary of named PETScMatrix and PETScVector objects
        to a file
        """
        assert self.simulation.ncpu == 1, 'Not supported in parallel'

        data = {}
        for key, value in kwargs.items():
            # Vectors
            if isinstance(value, dolfin.PETScVector):
                data[key] = ('dolfin.PETScVector', value.get_local())
            # Matrices
            elif isinstance(value, dolfin.PETScMatrix):
                rows = [0]
                cols = []
                values = []
                N, M = value.size(0), value.size(1)
                for irow in range(value.size(0)):
                    indices, row_values = value.getrow(irow)
                    rows.append(len(indices) + rows[-1])
                    cols.extend(indices)
                    values.extend(row_values)
                data[key] = ('dolfin.PETScMatrix', (N, M), rows, cols, values)

            else:
                raise ValueError('Cannot save object of type %r' % type(value))

        import pickle

        with open(file_name, 'wb') as out:
            pickle.dump(data, out, protocol=pickle.HIGHEST_PROTOCOL)
        self.simulation.log.info(
            'Saved LA objects to %r (%r)' % (file_name, kwargs.keys())
        )

    def load_la_objects(self, file_name):
        """
        Load a dictionary of named PETScMatrix and PETScVector objects
        from a file
        """
        assert self.simulation.ncpu == 1, 'Not supported in parallel'

        import pickle

        with open(file_name, 'rb') as inp:
            data = pickle.load(inp)

        ret = {}
        for key, value in data.items():
            value_type = value[0]
            # Vectors
            if value_type == 'dolfin.PETScVector':
                arr = value[1]
                dolf_vec = dolfin.PETScVector(dolfin.MPI.comm_world, arr.size)
                dolf_vec.set_local(arr)
                dolf_vec.apply('insert')
                ret[key] = dolf_vec
            # Matrices
            elif value_type == 'dolfin.PETScMatrix':
                shape, rows, cols, values = value[1:]
                from petsc4py import PETSc

                mat = PETSc.Mat().createAIJ(size=shape, csr=(rows, cols, values))
                mat.assemble()
                dolf_mat = dolfin.PETScMatrix(mat)
                ret[key] = dolf_mat

            else:
                raise ValueError('Cannot save object of type %r' % value_type)

        self.simulation.log.info(
            'Loaded LA objects from %r (%r)' % (file_name, data.keys())
        )
        return ret
