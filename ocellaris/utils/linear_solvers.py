import dolfin
import contextlib


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
    
    return make_linear_solver(solver_method, preconditioner, lu_method, params)


def make_linear_solver(solver_method, preconditioner=None, lu_method=None, parameters=None):
    """
    Create a Krylov or LU solver
    
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
    """
    if solver_method.lower() == 'lu':
        solver = dolfin.PETScLUSolver(lu_method)
    else:
        precon = dolfin.PETScPreconditioner(preconditioner)
        solver = dolfin.PETScKrylovSolver(solver_method, precon)
        solver.prec = precon # Keep from going out of scope
    
    for parameter_set in parameters:
        apply_settings(solver.parameters, parameter_set)
    
    solver.created_with_preconditioner = preconditioner
    solver.created_with_lu_method = lu_method
    solver.created_with_parameters = parameters
    
    return solver


def apply_settings(parameters, new_values):
    """
    This function does almost the same as::
    
        parameters.update(new_values)
        
    The difference is that subdictionaries are handled
    recursively and not replaced outright
    """
    for key, value in new_values.items():
        if isinstance(value, dict):
            apply_settings(parameters[key], value)
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
    for key, val in opts.iteritems():
        PETSc.Options().setValue(key, val)
    
    yield # run the code
    
    for key in opts.iterkeys():
        if key in orig_opts:
            PETSc.Options().setValue(key, orig_opts[key])
        else:
            PETSc.Options().delValue(key)


def condition_number(A):
    """
    Calculate the condition number of the given PetscMatrix A by use
    of SLEPSc solvers
    """
    from petsc4py import PETSc
    from slepc4py import SLEPc

    # Get the petc4py matrix
    PA = dolfin.as_backend_type(A).mat()

    # Calculate the largest and smallest singular value
    opts = {
        'svd_type': 'cross',
        'svd_eps_type': 'gd',
        'svd_converged_reason': 'ascii::ascii_info_detail',
        'svd_implicittranspose': True,
        'svd_tol': 1e-4,
        'svd_eps_refined': True,
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
            raise ValueError('Could not find the highest singular value')
        print 'Highest singular value:', sigma_1
        
        S.setWhichSingularTriplets(SLEPc.SVD.Which.SMALLEST)
        S.solve()
        if S.getConverged() == 1:
            sigma_n = S.getSingularTriplet(0)
        else:
            #raise ValueError('Could not find the lowest singular value')
            print S.getConvergedReason()
            sigma_n = 1
        print 'Lowest singular value:', sigma_n
        print PETSc.Options().getAll()
    print PETSc.Options().getAll()
    
    exit()
    return sigma_1/sigma_n
