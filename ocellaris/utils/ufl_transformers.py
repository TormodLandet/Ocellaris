import numpy
import dolfin
from ufl import as_vector, Form
from ufl.classes import FixedIndex, Indexed, ListTensor, MultiIndex, Zero
from ufl.algorithms import (expand_indices, expand_compounds,expand_derivatives,
                            compute_form_lhs, compute_form_rhs)
from ufl.corealg.multifunction import MultiFunction
from ufl.corealg.map_dag import map_expr_dag


#################################################
# Utility functions:


def split_form_into_matrix(full_form, Wv, Wu, Sv=None, Su=None, empty_cell_value=None, check_zeros=True):
    """
    Split a form into subforms which correspond to separate
    test and trial functions. Given a full form with multiple
    test and trial functions this will return a matrix of bi-
    linear forms and a vector of linear forms ordered in the
    same way as the ordering of the input test and trial
    functions.
    
    Given test functions (v, q) and trial functions (u, p)
    the return values are is a tuple consisting of the two
    lists:
    
    * Bilinear forms: [[A(u,v), B(p, v)], [C(u,q), D(p,q)]]
    * Linear forms: [E(v), F(q)]
    """
    N = Wv.num_sub_spaces()
    M = Wu.num_sub_spaces()
    form_matrix = numpy.zeros((N, M), dtype=numpy.object)
    form_vector = numpy.zeros(N, dtype=numpy.object)
    form_matrix[:] = form_vector[:] = empty_cell_value
    
    for i in range(N):
        # Process linear form
        f = FormPruner(i).prune(full_form, check_zeros)
        if f.integrals():
            form_vector[i] = f
        
        # Process bilinear form
        for j in range(M):
            f = FormPruner(i, j).prune(full_form, check_zeros)
            if f.integrals():
                form_matrix[i, j] = f
    
    return form_matrix, form_vector


def is_zero_ufl_expression(expr, return_val=False):
    """
    Is the given expression always identically zero or not
    Returns a boolean by default, but will return the actual
    evaluated expression value if return_val=True
    """
    # Reduce the complexity of the expression as much as possible
    expr = expand_derivatives(expr)
    expr = expand_compounds(expr)
    expr = expand_indices(expr)
    expr = IndexSimplificator().visit(expr)
    
    # Perform the zero-form estimation
    val = EstimateZeroForms().visit(expr)
    val = int(val)
    
    # val > 0 if the form is (likely) non-Zero and 0 if it is 
    # provably identically Zero()
    if return_val:
        return val
    else:
        return val == 0


#################################################
# Multifunction classes:


class FormPruner(MultiFunction):
    """
    Return a modified form where all arguments containing test
    or trial function with coupled function space indices which
    are not the specified indices (index_test & index_trial) are
    pruned from the UFL form expression tree
    
    You can use the "prune" method to create a pruned form
    """
    expr = MultiFunction.reuse_if_untouched
    
    def visit(self, expr):
        return map_expr_dag(self, expr)
    
    def __init__(self, index_test, index_trial=None):
        super(FormPruner, self).__init__()
        self._index_v = index_test
        self._index_u = index_trial
        self._cache = {}
    
    def prune(self, form, check_zeros=True):        
        # Get the parts of the form with the correct arity
        if self._index_u is None:
            form = compute_form_rhs(form)
        else:
            form = compute_form_lhs(form)
        
        integrals = []
        for integral in form.integrals():
            # Prune integrals that do not contain Arguments with
            # the chosen coupled function space indices 
            pruned = self.visit(integral.integrand())
            if not check_zeros or not is_zero_ufl_expression(pruned):
                integrals.append(integral.reconstruct(pruned))
        
        return Form(integrals)
    
    def argument(self, arg):
        """
        Argument is UFL lingo for test (num=0) and trial (num=1) functions
        """
        # Do not make several new args for the same original arg
        if arg in self._cache:
            return self._cache[arg]
        
        # Get some data about the argument and get our target index
        V = arg.function_space()
        N = V.num_sub_spaces()
        num = arg.number()
        idx_wanted = [self._index_v, self._index_u][num]
        
        new_arg = []
        for idx in range(N):
            # Construct non-coupled argument
            Vi = V.sub(idx).collapse()
            a = dolfin.function.argument.Argument(Vi, num, part=arg.part())
            indices = numpy.ndindex(a.ufl_shape)
            
            if idx == idx_wanted:
                # This is a wanted index
                new_arg.extend(a[I] for I in indices)
            else:
                # This index should be pruned
                new_arg.extend(Zero() for _ in indices)
        
        new_arg = as_vector(new_arg)
        self._cache[arg] = new_arg
        return new_arg


class IndexSimplificator(MultiFunction):
    """
    Simplifies indexing into ListTensors with fixed indices.
    from https://github.com/firedrakeproject/tsfc/compare/index-simplificator?expand=1#diff-a766247c71abcaca1251147d24ca9b63
    """
    expr = MultiFunction.reuse_if_untouched

    def visit(self, expr):
        return map_expr_dag(self, expr)

    def indexed(self, o, expr, multiindex):
        indices = list(multiindex)
        while indices and isinstance(expr, ListTensor) and isinstance(indices[0], FixedIndex):
            index = indices.pop(0)
            expr = expr.ufl_operands[int(index)]

        if indices:
            return Indexed(expr, MultiIndex(tuple(indices)))
        else:
            return expr


class EstimateZeroForms(MultiFunction):
    """
    Replace all non-Zero leaf nodes with 1 and then interpret the
    operator tree to calculate the scalar value of an expression
    in order to estimate if an UFL expression is allways identically
    zero or not.
    
    The value returned is the evaluated/interpreted expression. The 
    actual value of a non-zero return is not interesting in itself,
    but it indicates that the expression is potentially not identically
    zero when compiled by the form compiler.
    
    We aim to NEVER classify a non-zero form as Zero while detecting
    as many true zero forms as we can
    """
    def visit(self, expr):
        if isinstance(expr, (int, float)):
            return 1 if expr else 0
        return map_expr_dag(self, expr)

    # --- Terminal objects ---

    def zero(self, o):
        return 0

    def identity(self, o):
        raise NotImplementedError()

    def permutation_symbol(self, o):
        raise NotImplementedError()

    def facet_normal(self, o):
        assert len(o.ufl_shape) == 1
        return as_vector([1] * o.ufl_shape[0])

    def argument(self, o):
        shape = o.ufl_shape
        if len(shape) == 0:
            return 1
        elif len(shape) == 1:
            return as_vector([1]*shape[0])
        elif len(shape) == 2:
            return as_vector([as_vector([1]*shape[1]) for _ in range(shape[0])])
        else:
            raise NotImplementedError()
    constant = coefficient = scalar_value = argument

    def multi_index(self, o):
        return o # Handle further up

    def variable(self, o):
        raise NotImplementedError()

    # --- Non-terminal objects ---

    def index_sum(self, o, f, i):
        raise NotImplementedError()

    def sum(self, o, *ops):
        return sum(self.visit(o) for o in ops)

    def product(self, o, *ops):
        a = 1
        for oi in ops:
            a *= self.visit(oi)
        return a

    def division(self, o, a, b):
        return self.visit(a)

    def abs(self, o, a):
        return abs(self.visit(a))

    def transposed(self, o, a):
        raise NotImplementedError()

    def indexed(self, o, expr, multi_index):
        assert isinstance(multi_index, MultiIndex)
        indices = list(multi_index)
        
        if len(indices) == 1:
            assert isinstance(indices[0], FixedIndex)
            i = int(indices[0])
            return expr[i]
        
        elif len(indices) == 2:
            assert isinstance(indices[0], FixedIndex)
            assert isinstance(indices[1], FixedIndex)
            i = int(indices[0])
            j = int(indices[1])
            return expr[i][j]
        
        else:
            raise NotImplementedError()
    
    def variable_derivative(self, o, f, v):
        raise NotImplementedError()
    
    def coefficient_derivative(self, o, f, w, v):
        raise NotImplementedError()
    
    def grad(self, o, expr):
        shape = o.ufl_shape
        if len(shape) == 1:
            return as_vector([expr] * shape[0])
        elif len(shape) == 2:
            return as_vector([[expr[i]] * shape[1] for i in range(shape[0])])
        else:
            raise NotImplementedError()
    
    def div(self, o, expr):
        shape = o.ufl_shape
        if len(shape) == 0:
            return numpy.max(expr)
        else:
            raise NotImplementedError()
    
    def nabla_grad(self, o, f):
        raise NotImplementedError()

    def nabla_div(self, o, f):
        raise NotImplementedError()

    def curl(self, o, f):
        raise NotImplementedError()

    # Functions of one scalar argument that are zero for zero arguments
    def sqrt(self, o, f):
        return self.visit(f)
    sin = tan = errf = sqrt
    sinh = tanh = asin = atan = sqrt
    
    # Functions of one scalar argument that are non-zero for zero arguments
    def cos(self, o, f):
        return 1
    ln = exp = cosh = acos = cos
    
    # Functions of two scalar arguments
    def atan2(self, o, f1, f2):
        return 1
    bessel_j = bessel_y = bessel_i = bessel_K = atan2
    
    def power(self, o, a, b):
        return self.visit(a)

    def outer(self, o, a, b):
        raise NotImplementedError()

    def inner(self, o, a, b):
        raise NotImplementedError()

    def dot(self, o, a, b):
        return numpy.dot(a, b)
    
    def cross(self, o, a, b):
        raise NotImplementedError()

    def trace(self, o, A):
        raise NotImplementedError()

    def determinant(self, o, A):
        raise NotImplementedError()

    def inverse(self, o, A):
        raise NotImplementedError()

    def deviatoric(self, o, A):
        raise NotImplementedError()

    def cofactor(self, o, A):
        raise NotImplementedError()

    def skew(self, o, A):
        raise NotImplementedError()

    def sym(self, o, A):
        raise NotImplementedError()

    def list_tensor(self, o):
        shape = o.ufl_shape
        if len(shape) == 1:
            return as_vector([self.visit(op) for op in o.ufl_operands])
        elif len(shape) == 2:
            return as_vector([as_vector([self.visit(op) for op in row.ufl_operands])
                              for row in o.ufl_operands])
        else:
            raise NotImplementedError()

    def component_tensor(self, o, *ops):
        raise NotImplementedError()

    def positive_restricted(self, o, f):
        return self.visit(f)

    negative_restricted = cell_avg = positive_restricted

    # The value of a condition is not important, 
    # we will assume both true and false anyway
    def eq(self, o, a, b):
        return 1
    def not_condition(self, o, a):
        raise 1 
    ne = le = ge = lt = and_condition = or_condition = eq
    def conditional(self, o, c, t, f):
        if o.ufl_shape == ():
            return max(t, f)
        else:
            raise NotImplementedError()

    def min_value(self, o, a, b):
        a = self.visit(a)
        b = self.visit(b)
        return max(a, b) # not a typo, should be max!

    max_value = min_value

    def expr(self, o):
        raise NotImplementedError("Missing handler for type %s" % str(type(o)))
