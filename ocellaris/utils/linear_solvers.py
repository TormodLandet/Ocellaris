import numpy
import dolfin
import contextlib
from petsc4py import PETSc
from ocellaris.cpp import load_module


def linear_solver_from_input(simulation, path, default_solver, default_preconditioner,
                             default_lu_method, default_parameters):
    """
    From specifications in the input at the given path create a linear solver
    
    The path (e.g "solver/u") must point to a dictionary in the input file that
    can contain optional fields specifying the solver.
    
    Example::
    
        solver:
            u:
                solver: gmres
                preconditioner: additive_schwarz
            coupled:
                solver: lu
                lu_method: mumps
                parameters:
                    same_nonzero_pattern: True
    
    The default values are used if the keys are not found in the input
    """
    # Get values from input dictionary
    solver_method = simulation.input.get_value('%s/solver' % path, default_solver, 'string')
    preconditioner = simulation.input.get_value('%s/preconditioner' % path, default_preconditioner, 'string')
    lu_method = simulation.input.get_value('%s/lu_method' % path, default_lu_method, 'string')
    solver_parameters = simulation.input.get_value('%s/parameters' % path, {}, 'dict(string:any)')
    params = [default_parameters, solver_parameters]
    
    simulation.log.info('    Creating linear equation solver from input "%s"' % path)
    simulation.log.info('        Method:         %s' % solver_method)
    simulation.log.info('        Preconditioner: %s' % preconditioner)
    simulation.log.info('        LU-method:      %s' % lu_method)
    
    return LinearSolverWrapper(solver_method, preconditioner, lu_method, params)


class LinearSolverWrapper(object):
    def __init__(self, solver_method, preconditioner=None, lu_method=None, parameters=None): 
        """
        Wrap a Krylov or LU solver
        
        You must either specify solver_method = 'lu' and give the name
        of the solver, e.g lu_solver='mumps' or give a valid Krylov
        solver name, eg. solver_method='minres' and give the name of a
        preconditioner, eg. preconditioner_name='hypre_amg'.
        
        The parameters argument is a *list* of dictionaries which are
        to be used as parameters to the Krylov solver. Settings in the
        first dictionary in this list will be (potentially) overwritten
        by settings in later dictionaries. The use case is to provide
        sane defaults as well as allow the user to override the defaults
        in the input file
        
        The reason for this wrapper is to provide easy querying of 
        iterative/direct and not crash when set_reuse_preconditioner is 
        run before the first solve. This simplifies usage
        """
        self.solver_method = solver_method
        self.preconditioner = preconditioner
        self.lu_method = lu_method
        self.input_parameters = parameters
        
        self.is_first_solve = True
        self.is_iterative = False
        self.is_direct = False
    
        if solver_method.lower() == 'lu':
            solver = dolfin.PETScLUSolver(lu_method)
            self.is_direct = True
        else:
            precon = dolfin.PETScPreconditioner(preconditioner)
            solver = dolfin.PETScKrylovSolver(solver_method, precon)
            self._pre_obj = precon # Keep from going out of scope
            self.is_iterative = True
        
        for parameter_set in parameters:
            apply_settings(solver_method, solver.parameters, parameter_set)
        
        self._solver = solver
    
    def solve(self, *argv, **kwargs):
        ret = self._solver.solve(*argv, **kwargs)
        self.is_first_solve = False
        return ret
    
    @property
    def parameters(self):
        return self._solver.parameters 
    
    def set_reuse_preconditioner(self, *argv, **kwargs):
        if self.is_iterative and self.is_first_solve:
            return  # Nov 2016: this segfaults if running before the first solve
        else:
            return self._solver.set_reuse_preconditioner(*argv, **kwargs)
    
    def __repr__(self):
        return ('<LinearSolverWrapper iterative=%r ' % self.is_iterative + 
                                     'direct=%r ' % self.is_direct +
                                     'method=%r ' % self.solver_method +
                                     'preconditioner=%r ' % self.preconditioner +
                                     'LU-method=%r ' % self.lu_method +
                                     'parameters=%r>' % self.input_parameters)


def apply_settings(solver_method, parameters, new_values):
    """
    This function does almost the same as::
    
        parameters.update(new_values)
        
    The difference is that subdictionaries are handled
    recursively and not replaced outright
    """
    skip = set()
    if solver_method == 'lu':
        skip.update(['nonzero_initial_guess',
                     'relative_tolerance',
                     'absolute_tolerance'])
    
    for key, value in new_values.items():
        if key in skip:
            continue
        elif isinstance(value, dict):
            apply_settings(solver_method, parameters[key], value)
        else:
            parameters[key] = value


@contextlib.contextmanager
def petsc_options(opts):
    """
    A context manager to set PETSc options for a limited amount of code.
    The parameter opts is a dictionary of PETSc/SLEPc options 
    """
    from petsc4py import PETSc
    orig_opts = PETSc.Options().getAll()
    for key, val in opts.items():
        PETSc.Options().setValue(key, val)
    
    yield # run the code
    
    for key in opts.keys():
        if key in orig_opts:
            PETSc.Options().setValue(key, orig_opts[key])
        else:
            PETSc.Options().delValue(key)


def create_block_matrix(V, blocks):
    """
    Create a sparse matrix to hold dense blocks that are larger than
    the normal DG block diagonal mass matrices (super-cell dense blocks)
    
    The argument ``blocks`` should be a list of lists/arrays containing
    the dofs in each block. The dofs are assumed to be the same for
    both rows and columns.
    """
    comm = V.mesh().mpi_comm()
    dm = V.dofmap()
    im = dm.index_map()
    
    # Create a tensor layout for the matrix
    ROW_MAJOR = 0
    tl = dolfin.TensorLayout(comm, ROW_MAJOR, dolfin.TensorLayout.Sparsity.SPARSE)
    tl.init([im, im], dolfin.TensorLayout.Ghosts.GHOSTED)
    
    # Setup the tensor layout's sparsity pattern
    sp = tl.sparsity_pattern()
    sp.init([im, im])
    entries = None
    for block in blocks:
        N = len(block)
        if entries is None or entries.shape[1] != N:
            entries = numpy.empty((2, N), dtype=numpy.intc)
        entries[0,:] = block
        entries[1,:] = entries[0,:]
        sp.insert_local(entries)
    sp.apply()
    
    # Create a matrix with the newly created tensor layout
    A = dolfin.PETScMatrix(comm)
    A.init(tl)
    
    return A


def matmul(A, B, out=None, use_cpp=False):
    """
    A B (and potentially out) must be PETScMatrix
    The matrix out must be the result of a prior matmul
    call with the same sparsity patterns in A and B
    """
    assert A is not None and B is not None
    
    if use_cpp:
        # Implemented in C++ due to unfinished petsc4py support in 
        # dolfin's new pybind11 version (as of 2017-09-25)
        mod = load_module('petsc_utils')
        if out is None:
            C = mod.matmul(A, B)
        else:
            C = out
            mod.matmul_reuse(A, B, C)
    
    else:
        A = A.mat()
        B = B.mat()
        if out is not None:
            A.matMult(B, out.mat())
            C = out
        else:
            Cmat = A.matMult(B)
            C = dolfin.PETScMatrix(Cmat)
            C.apply('insert')
    
    return C


def condition_number(A, method='simplified'):
    """
    Estimate the condition number of the matrix A
    """
    if method == 'simplified':
        # Calculate max(abs(A))/min(abs(A))
        amin, amax = 1e10, -1e10
        for irow in range(A.size(0)):
            _indices, values = A.getrow(irow)
            aa = abs(values)
            amax = max(amax, aa.max())
            aa[aa==0] = amax
            amin = min(amin, aa.min())
        amin = dolfin.MPI.min(dolfin.mpi_comm_world(), float(amin))
        amax = dolfin.MPI.max(dolfin.mpi_comm_world(), float(amax))
        return amax/amin
    
    elif method == 'numpy':
        from numpy.linalg import cond
        A = mat_to_scipy_csr(A).todense()
        return cond(A)
    
    elif method == 'SLEPc':
        from slepc4py import SLEPc
        
        # Get the petc4py matrix
        PA = dolfin.as_backend_type(A).mat()
        
        # Calculate the largest and smallest singular value
        opts = {
            'svd_type': 'cross',
            'svd_eps_type': 'gd',
            #'help': 'svd_type'
        }
        with petsc_options(opts):
            S = SLEPc.SVD()
            S.create()
            S.setOperator(PA)
            S.setFromOptions()
            S.setDimensions(1, PETSc.DEFAULT, PETSc.DEFAULT)
            S.setWhichSingularTriplets(SLEPc.SVD.Which.LARGEST)
            S.solve()
            if S.getConverged() == 1:
                sigma_1 = S.getSingularTriplet(0)
            else:
                raise ValueError('Could not find the highest singular value (%d)'
                                 % S.getConvergedReason())
            print('Highest singular value:', sigma_1)
            
            S.setWhichSingularTriplets(SLEPc.SVD.Which.SMALLEST)
            S.solve()
            if S.getConverged() == 1:
                sigma_n = S.getSingularTriplet(0)
            else:
                raise ValueError('Could not find the lowest singular value (%d)'
                                 % S.getConvergedReason())
            print('Lowest singular value:', sigma_n)
            print(PETSc.Options().getAll())
        print(PETSc.Options().getAll())
        
        return sigma_1/sigma_n


def mat_to_scipy_csr(dolfin_matrix):
    """
    Convert any dolfin.Matrix to csr matrix in scipy.
    Based on code by Miroslav Kuchta
    """
    assert dolfin.MPI.size(dolfin.mpi_comm_world()) == 1, 'mat_to_csr assumes single process'
    import scipy.sparse
    import numpy
    
    rows = [0]
    cols = []
    values = []
    for irow in range(dolfin_matrix.size(0)):
        indices, values_ = dolfin_matrix.getrow(irow)
        rows.append(len(indices)+rows[-1])
        cols.extend(indices)
        values.extend(values_)

    shape = dolfin_matrix.size(0), dolfin_matrix.size(1)
        
    return scipy.sparse.csr_matrix((numpy.array(values, dtype='float'),
                                    numpy.array(cols, dtype='int'),
                                    numpy.array(rows, dtype='int')),
                                    shape)
