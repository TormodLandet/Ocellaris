import dolfin
import numpy
from ocellaris.utils import timeit

@timeit
def velocity_error_norm(u1, u2, Vu_highp, norm='L2'):
    """
    Compute the difference in the given norm between the 
    two velocity fields u1 and u2
    """
    ndim = Vu_highp.cell().topological_dimension()
    
    if norm == 'L2':
        # Interpolate velocities to higher polynomial space
        u1_vec = []
        u2_vec = []
        for d in range(ndim):
            u1_vec.append(dolfin.interpolate(u1[d], Vu_highp))
            u2_vec.append(dolfin.interpolate(u2[d], Vu_highp))
        
        # Recreate vector valued function
        u1_highp = dolfin.as_vector(u1_vec)
        u2_highp = dolfin.as_vector(u2_vec)
        
        # Calculate the L^2 norm
        res = dolfin.assemble((u1_highp - u2_highp)**2*dolfin.dx)**0.5
    
    elif norm == 'Linf':
        # Falculate L^inf norm
        res = 0
        for d in range(ndim):
            a1 = u1[d].vector().array()
            a2 = u2[d].vector().array()
            err = numpy.abs(a1 - a2).max()  
            res = max(res, err)
    
    return res
