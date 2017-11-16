"""
Convert to and from a Taylor basis

This code follows the definition of cell average Taylor DG polynomial expansion
as described in Dmitri Kuzmin (2010) "A vertex-based hierarchial slope limiter
for p-adaptive discontinuous Galerkin methods"
"""
import numpy
from dolfin import cells
from ocellaris.utils import ocellaris_error


CACHE = {}


def lagrange_to_taylor(u, t):
    """
    Convert a Lagrange function space function into a DG Taylor
    function space function. The results are stored in t which
    should be a function space with an appropriate number of dofs
    per cell (t in DG2 if u in DG2 etc)
    
    Note: Most possible dolfin operations on the function t will
    be wrong since dolfin does not support Taylor basis DG elements.
    The t function should hence be read and modified by specially
    crafted code only! 
    """
    V = u.function_space()
    mesh = V.mesh()
    degree = V.ufl_element().degree()
    ndim = mesh.geometry().dim()
    
    if 'mesh_hash' not in CACHE or CACHE['mesh_hash'] != mesh.hash():
        CACHE.clear()
        CACHE['mesh'] = mesh
        CACHE['mesh_hash'] = mesh.hash()
    
    # Get cached conversion matrices
    key = ('lagrange_to_taylor_matrices', degree)
    if key not in CACHE:
        if degree == 1 and ndim == 2:
            CACHE[key] = DG1_to_taylor_matrix_2D(V)
        elif degree == 1 and ndim == 3:
            CACHE[key] = DG1_to_taylor_matrix_3D(V)
        elif degree == 2 and ndim == 2:
            CACHE[key] = DG2_to_taylor_matrix_2D(V)
        elif degree == 2 and ndim == 3:
            CACHE[key] = DG2_to_taylor_matrix_3D(V)
        else:
            ocellaris_error('DG Lagrange to DG Taylor converter error',
                            'Polynomial degree %d not supported' % degree)
    lagrange_to_taylor_matrices = CACHE[key]
    Ncells = lagrange_to_taylor_matrices.shape[0]
    
    # Get cached cell dofs
    key2 = ('cell_dofs', degree)
    if key2 not in CACHE:
        dm = V.dofmap()
        cell_dofs = [dm.cell_dofs(i) for i in range(Ncells)]
        CACHE[key2] = numpy.array(cell_dofs, int)
    cell_dofs = CACHE[key2]
    
    # Apply the conversion matrices for all cells by use of the stacked dot
    # behaviour of matmul (slightly faster than einsum 'ijk,ik->ij')
    all_vals_lagrange = u.vector().get_local()
    lagrange_vectors = all_vals_lagrange.take(cell_dofs)
    res = numpy.matmul(lagrange_to_taylor_matrices, lagrange_vectors[:,:,None]).squeeze()
    
    # Put the results into the right indices in the Taylor function's vector
    all_vals_taylor = numpy.zeros_like(all_vals_lagrange)
    all_vals_taylor[cell_dofs] = res
    t.vector().set_local(all_vals_taylor)
    t.vector().apply('insert')


def taylor_to_lagrange(t, u):
    """
    Convert a DG Taylor function space function into a Lagrange
    function space function. The results are stored in u which
    should be a function space with an appropriate number of dofs
    per cell (u in DG2 if t in DG2 etc)
    """
    V = u.function_space()
    mesh = V.mesh()
    degree = V.ufl_element().degree()
    ndim = mesh.geometry().dim()
    
    if 'mesh_hash' not in CACHE or CACHE['mesh_hash'] != mesh.hash():
        CACHE.clear()
        CACHE['mesh'] = mesh
        CACHE['mesh_hash'] = mesh.hash()
    
    # Get cached conversion matrices
    key = ('taylor_to_lagrange_matrices', degree)
    if key not in CACHE:
        if degree == 1 and ndim == 2:
            CACHE[key] = taylor_to_DG1_matrix_2D(V)
        elif degree == 1 and ndim == 3:
            CACHE[key] = taylor_to_DG1_matrix_3D(V)
        elif degree == 2 and ndim == 2:
            CACHE[key] = taylor_to_DG2_matrix_2D(V)
        elif degree == 2 and ndim == 3:
            CACHE[key] = taylor_to_DG2_matrix_3D(V)
        else:
            ocellaris_error('DG Taylor to DG Lagrange converter error',
                            'Polynomial degree %d not supported' % degree)
    taylor_to_lagrange_matrices = CACHE[key]
    Ncells = taylor_to_lagrange_matrices.shape[0]
    
    # Get cached cell dofs
    key2 = ('cell_dofs', degree)
    if key2 not in CACHE:
        dm = V.dofmap()
        cell_dofs = [dm.cell_dofs(i) for i in range(Ncells)]
        CACHE[key2] = numpy.array(cell_dofs, int)
    cell_dofs = CACHE[key2]
    
    # Apply the conversion matrices for all cells by use of the stacked dot
    # behaviour of matmul (slightly faster than einsum 'ijk,ik->ij')
    all_vals_taylor = t.vector().get_local()
    taylor_vectors = all_vals_taylor.take(cell_dofs)
    res = numpy.matmul(taylor_to_lagrange_matrices, taylor_vectors[:,:,None]).squeeze()
    
    # Put the results into the right indices in the Taylor function's vector
    all_vals_lagrange = numpy.zeros_like(all_vals_taylor)
    all_vals_lagrange[cell_dofs] = res
    u.vector().set_local(all_vals_lagrange)
    u.vector().apply('insert')


##############################################################################


def DG1_to_taylor_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG1 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()
    
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    
    A = numpy.zeros((num_cells_owned, 3, 3), float)
    for cell in cells(mesh):
        icell = cell.index()
        
        verts = cell.entities(0)
        x = numpy.zeros((3, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
                
        ###############################
        # From sympy code gen
        
        ((x1, y1), (x2, y2), (x3, y3)) = x
        D = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
        
        # Value at xc (also the cell average value)
        A[icell,0,0] = 1/3
        A[icell,0,1] = 1/3
        A[icell,0,2] = 1/3
        
        # d/dx
        A[icell,1,0] = (y2 - y3)/D
        A[icell,1,1] = (-y1 + y3)/D
        A[icell,1,2] = (y1 - y2)/D
        
        # d/dy
        A[icell,2,0] = (-x2 + x3)/D
        A[icell,2,1] = (x1 - x3)/D
        A[icell,2,2] = (-x1 + x2)/D

    return A


def DG1_to_taylor_matrix_3D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG1 function space V in 3D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()
    
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    
    A = numpy.zeros((num_cells_owned, 4, 4), float)
    for cell in cells(mesh):
        icell = cell.index()
        
        verts = cell.entities(0)
        x = numpy.zeros((4, 3), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = vertices[verts[3]]
                
        ###############################
        # From sympy code gen
        
        ((x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4)) = x
        
        V = abs((-x1 + x2)*((-y1 + y3)*(-z1 + z4) - (-y1 + y4)*(-z1 + z3))/6 +
                (-y1 + y2)*(-(-x1 + x3)*(-z1 + z4) + (-x1 + x4)*(-z1 + z3))/6 +
                (-z1 + z2)*((-x1 + x3)*(-y1 + y4) - (-x1 + x4)*(-y1 + y3))/6)

        # Value at xc (also the cell average value)
        A[icell,0,0] = 1/4
        A[icell,0,1] = 1/4
        A[icell,0,2] = 1/4
        A[icell,0,3] = 1/4
        
        # Value at xc (also the cell average value)
        A[icell,0,0] = 1/4
        A[icell,0,1] = 1/4
        A[icell,0,2] = 1/4
        A[icell,0,3] = 1/4
        
        def sign(a):
            return 1 if a >= 0 else -1
        
        # d/dx
        A[icell,1,0] = ((y2 - y3)*(z2 - z4) - (y2 - y4)*(z2 - z3))*sign(((-x2 + x3)*(-y2 + y4) - (-x2 + x4)*(-y2 + y3))*(z1 - 3*z2 + z3 + z4) - ((x2 - x3)*(z2 - z4) - (x2 - x4)*(z2 - z3))*(y1 - 3*y2 + y3 + y4) + ((-y2 + y3)*(-z2 + z4) - (-y2 + y4)*(-z2 + z3))*(x1 - 3*x2 + x3 + x4))/(6*V)
        A[icell,1,1] = ((y1 - y3)*(z1 - z4) - (y1 - y4)*(z1 - z3))*sign(((x1 - x3)*(y1 - y4) - (x1 - x4)*(y1 - y3))*(-3*z1 + z2 + z3 + z4) - ((x1 - x3)*(z1 - z4) - (x1 - x4)*(z1 - z3))*(-3*y1 + y2 + y3 + y4) + ((y1 - y3)*(z1 - z4) - (y1 - y4)*(z1 - z3))*(-3*x1 + x2 + x3 + x4))/(6*V)
        A[icell,1,2] = ((y1 - y2)*(z1 - z4) - (y1 - y4)*(z1 - z2))*sign(((-x1 + x2)*(-y1 + y4) - (-x1 + x4)*(-y1 + y2))*(-3*z1 + z2 + z3 + z4) - ((x1 - x2)*(z1 - z4) - (x1 - x4)*(z1 - z2))*(-3*y1 + y2 + y3 + y4) + ((-y1 + y2)*(-z1 + z4) - (-y1 + y4)*(-z1 + z2))*(-3*x1 + x2 + x3 + x4))/(6*V)
        A[icell,1,3] = ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*sign(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*(-3*z1 + z2 + z3 + z4) - ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*(-3*y1 + y2 + y3 + y4) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*(-3*x1 + x2 + x3 + x4))/(6*V)
        
        # d/dy
        A[icell,2,0] = -((x2 - x3)*(z2 - z4) - (x2 - x4)*(z2 - z3))*sign(((x2 - x3)*(y2 - y4) - (x2 - x4)*(y2 - y3))*(z1 - 3*z2 + z3 + z4) - ((x2 - x3)*(z2 - z4) - (x2 - x4)*(z2 - z3))*(y1 - 3*y2 + y3 + y4) + ((y2 - y3)*(z2 - z4) - (y2 - y4)*(z2 - z3))*(x1 - 3*x2 + x3 + x4))/(6*V)
        A[icell,2,1] = -((x1 - x3)*(z1 - z4) - (x1 - x4)*(z1 - z3))*sign(((x1 - x3)*(y1 - y4) - (x1 - x4)*(y1 - y3))*(-3*z1 + z2 + z3 + z4) - ((x1 - x3)*(z1 - z4) - (x1 - x4)*(z1 - z3))*(-3*y1 + y2 + y3 + y4) + ((y1 - y3)*(z1 - z4) - (y1 - y4)*(z1 - z3))*(-3*x1 + x2 + x3 + x4))/(6*V)
        A[icell,2,2] = -((x1 - x2)*(z1 - z4) - (x1 - x4)*(z1 - z2))*sign(((x1 - x2)*(y1 - y4) - (x1 - x4)*(y1 - y2))*(-3*z1 + z2 + z3 + z4) - ((x1 - x2)*(z1 - z4) - (x1 - x4)*(z1 - z2))*(-3*y1 + y2 + y3 + y4) + ((y1 - y2)*(z1 - z4) - (y1 - y4)*(z1 - z2))*(-3*x1 + x2 + x3 + x4))/(6*V)
        A[icell,2,3] = -((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*sign(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*(-3*z1 + z2 + z3 + z4) - ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*(-3*y1 + y2 + y3 + y4) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*(-3*x1 + x2 + x3 + x4))/(6*V)
        
        # d/dz
        A[icell,3,0] = ((x2 - x3)*(y2 - y4) - (x2 - x4)*(y2 - y3))*sign(((-x2 + x3)*(-y2 + y4) - (-x2 + x4)*(-y2 + y3))*(z1 - 3*z2 + z3 + z4) - ((x2 - x3)*(z2 - z4) - (x2 - x4)*(z2 - z3))*(y1 - 3*y2 + y3 + y4) + ((-y2 + y3)*(-z2 + z4) - (-y2 + y4)*(-z2 + z3))*(x1 - 3*x2 + x3 + x4))/(6*V)
        A[icell,3,1] = ((x1 - x3)*(y1 - y4) - (x1 - x4)*(y1 - y3))*sign(((x1 - x3)*(y1 - y4) - (x1 - x4)*(y1 - y3))*(-3*z1 + z2 + z3 + z4) - ((x1 - x3)*(z1 - z4) - (x1 - x4)*(z1 - z3))*(-3*y1 + y2 + y3 + y4) + ((y1 - y3)*(z1 - z4) - (y1 - y4)*(z1 - z3))*(-3*x1 + x2 + x3 + x4))/(6*V)
        A[icell,3,2] = ((x1 - x2)*(y1 - y4) - (x1 - x4)*(y1 - y2))*sign(((-x1 + x2)*(-y1 + y4) - (-x1 + x4)*(-y1 + y2))*(-3*z1 + z2 + z3 + z4) - ((x1 - x2)*(z1 - z4) - (x1 - x4)*(z1 - z2))*(-3*y1 + y2 + y3 + y4) + ((-y1 + y2)*(-z1 + z4) - (-y1 + y4)*(-z1 + z2))*(-3*x1 + x2 + x3 + x4))/(6*V)
        A[icell,3,3] = ((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*sign(((x1 - x2)*(y1 - y3) - (x1 - x3)*(y1 - y2))*(-3*z1 + z2 + z3 + z4) - ((x1 - x2)*(z1 - z3) - (x1 - x3)*(z1 - z2))*(-3*y1 + y2 + y3 + y4) + ((y1 - y2)*(z1 - z3) - (y1 - y3)*(z1 - z2))*(-3*x1 + x2 + x3 + x4))/(6*V)
    
    return A

def DG2_to_taylor_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Lagrange cell dofs return a vector of Taylor cell dofs.
    This implementation handles DG2 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()
    
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    
    A = numpy.zeros((num_cells_owned, 6, 6), float)
    for cell in cells(mesh):
        icell = cell.index()
        
        verts = cell.entities(0)
        x = numpy.zeros((3, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
                
        ###############################
        # From sympy code gen
        
        ((x1, y1), (x2, y2), (x3, y3)) = x
        D = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
        
        # MODIFIED to give the average value, not the center value!
        A[icell,0,3] = 1/3
        A[icell,0,4] = 1/3
        A[icell,0,5] = 1/3
        
        # d/dx
        A[icell,1,0] = (y2 - y3)/(3*D)
        A[icell,1,1] = (-y1 + y3)/(3*D)
        A[icell,1,2] = (y1 - y2)/(3*D)
        A[icell,1,3] = 4*(-y2 + y3)/(3*D)
        A[icell,1,4] = 4*(y1 - y3)/(3*D)
        A[icell,1,5] = 4*(-y1 + y2)/(3*D)
        
        # d/dy
        A[icell,2,0] = (-x2 + x3)/(3*D)
        A[icell,2,1] = (x1 - x3)/(3*D)
        A[icell,2,2] = (-x1 + x2)/(3*D)
        A[icell,2,3] = 4*(x2 - x3)/(3*D)
        A[icell,2,4] = 4*(-x1 + x3)/(3*D)
        A[icell,2,5] = 4*(x1 - x2)/(3*D)
        
        # d/dx^2
        A[icell,3,0] = 4*(y2 - y3)**2/D**2
        A[icell,3,1] = 4*(y1 - y3)**2/D**2
        A[icell,3,2] = 4*(y1 - y2)**2/D**2
        A[icell,3,3] = -8*(y1 - y2)*(y1 - y3)/D**2
        A[icell,3,4] = 8*(y1 - y2)*(y2 - y3)/D**2
        A[icell,3,5] = -8*(y1 - y3)*(y2 - y3)/D**2
        
        # d/dy^2
        A[icell,4,0] = 4*(x2 - x3)**2/D**2
        A[icell,4,1] = 4*(x1 - x3)**2/D**2
        A[icell,4,2] = 4*(x1 - x2)**2/D**2
        A[icell,4,3] = -8*(x1 - x2)*(x1 - x3)/D**2
        A[icell,4,4] = 8*(x1 - x2)*(x2 - x3)/D**2
        A[icell,4,5] = -8*(x1 - x3)*(x2 - x3)/D**2
        
        # d/dx*dy
        A[icell,5,0] = -4*(x2 - x3)*(y2 - y3)/D**2
        A[icell,5,1] = -4*(x1 - x3)*(y1 - y3)/D**2
        A[icell,5,2] = -4*(x1 - x2)*(y1 - y2)/D**2
        A[icell,5,3] = 4*((x1 - x2)*(y1 - y3) + (x1 - x3)*(y1 - y2))/D**2
        A[icell,5,4] = -(4*(x1 - x2)*(y2 - y3) + 4*(x2 - x3)*(y1 - y2))/D**2
        A[icell,5,5] = 4*((x1 - x3)*(y2 - y3) + (x2 - x3)*(y1 - y3))/D**2

    return A


def taylor_to_DG1_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG1 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()
    
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    
    x = numpy.zeros((3, 2), float)
    A = numpy.zeros((num_cells_owned, 3, 3), float)
    for cell in cells(mesh):
        icell = cell.index()
        
        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        xc = (x[0] + x[1] +  x[2])/3
        
        for i in range(3):
            dx, dy = x[i,0] - xc[0], x[i,1] - xc[1]
            A[icell,i,0] = 1
            A[icell,i,1] = dx
            A[icell,i,2] = dy
    
    return A


def taylor_to_DG1_matrix_3D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG1 function space V in 3D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()
    
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    
    x = numpy.zeros((4, 3), float)
    A = numpy.zeros((num_cells_owned, 4, 4), float)
    for cell in cells(mesh):
        icell = cell.index()
        
        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = vertices[verts[3]]
        xc = (x[0] + x[1] +  x[2] +  x[3])/4
        
        for i in range(4):
            dx, dy, dz = x[i,0] - xc[0], x[i,1] - xc[1], x[i,2] - xc[2]
            A[icell,i,0] = 1
            A[icell,i,1] = dx
            A[icell,i,2] = dy
            A[icell,i,3] = dz
    
    return A


def taylor_to_DG2_matrix_2D(V):
    """
    Create the per cell matrices that when matrix multiplied with the
    Taylor cell dofs return a vector of Lagrange cell dofs.
    This implementation handles DG2 function space V in 2D
    """
    mesh = V.mesh()
    vertices = mesh.coordinates()
    
    tdim = mesh.topology().dim()
    num_cells_owned = mesh.topology().ghost_offset(tdim)
    
    x = numpy.zeros((6, 2), float)
    A = numpy.zeros((num_cells_owned, 6, 6), float)
    for cell in cells(mesh):
        icell = cell.index()
        
        verts = cell.entities(0)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = (x[1] + x[2])/2
        x[4] = (x[0] + x[2])/2
        x[5] = (x[0] + x[1])/2
        xc = (x[0] + x[1] +  x[2])/3
        
        # Code generated by the sympy code included below
        ((x1, y1), (x2, y2), (x3, y3)) = x[:3]
        bar_xx = x1**2/36 - x1*x2/36 - x1*x3/36 + x2**2/36 - x2*x3/36 + x3**2/36
        bar_yy = y1**2/36 - y1*y2/36 - y1*y3/36 + y2**2/36 - y2*y3/36 + y3**2/36
        bar_xy = (x1*y1/18 - x1*y2/36 - x1*y3/36 - x2*y1/36 + x2*y2/18 - x2*y3/36
                  - x3*y1/36 - x3*y2/36 + x3*y3/18)
        
        for i in range(6):
            dx, dy = x[i,0] - xc[0], x[i,1] - xc[1]
            A[icell,i,0] = 1
            A[icell,i,1] = dx
            A[icell,i,2] = dy
            A[icell,i,3] = dx**2/2 - bar_xx
            A[icell,i,4] = dy**2/2 - bar_yy
            A[icell,i,5] = dx*dy  - bar_xy
    
    return A


#########################################################################################
# Code production with SymPy
#
#   The code below produces the coefficients and expressions used above
#
#   This code is not used during normal operations, it is for development
#   only and is just kept here to help us remember that it exists and
#   does not have to be reinvented every time

def produce_code_with_sympy_DG1_2D():
    """
    Use sympy to derive the coefficients needed to convert from DG1 polynomials
    on a triangle to a Taylor basis by expressing the derivatives at the centre
    point of the triangle
    (Only used for generating the code above, not used while running Ocellaris)
    """
    from sympy import Symbol, symbols, S
    x1, x2, x3, y1, y2, y3 = symbols('x1 x2 x3 y1 y2 y3')
    x, y, D = symbols('x y D')
    
    # Barycentric coordinates from cartesian corner coordinates
    L1f = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))/D
    L2f = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))/D
    L3f = 1 - L1f - L2f
    L1, L2, L3 = Symbol('L1')(x, y), Symbol('L2')(x, y), Symbol('L3')(x, y)
    
    replacements = [(L1.diff(x),  L1f.diff(x)),
                    (L2.diff(x),  L2f.diff(x)),
                    (L3.diff(x),  L3f.diff(x)),
                    (L1.diff(y),  L1f.diff(y)),
                    (L2.diff(y),  L2f.diff(y)),
                    (L3.diff(y),  L3f.diff(y)),
                    (L1, 1/S(3)),
                    (L2, 1/S(3)),
                    (L3, 1/S(3))]
    
    # Barycentric quadratic shape functions on a triangle
    phi = [L1, L2, L3]
    
    def print_code(name, func, index):
        code = ['# %s' % name]
        for i in range(3):
            v = func(phi[i])
            v = v.subs(replacements)
            v = v.simplify()
            code.append('A[icell,%d,%d] = %s' % (index, i, v))
        code.append('')
        print('\n'.join(code))
    
    print_code('Value at xc (also the cell average value)', lambda q: q, 0)
    print_code('d/dx', lambda q: q.diff(x), 1)
    print_code('d/dy', lambda q: q.diff(y), 2)


def produce_code_with_sympy_DG1_3D(degree=1):
    """
    Use sympy to derive the coefficients needed to convert from DG1 polynomials
    on a tetrahedron to a Taylor basis by expressing the derivatives at the centre
    point of the triangle
    (Only used for generating the code above, not used while running Ocellaris)
    """
    from sympy import Symbol, symbols, S
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4', real=True)
    y1, y2, y3, y4 = symbols('y1 y2 y3 y4', real=True)
    z1, z2, z3, z4 = symbols('z1 z2 z3 z4', real=True)
    x, y, z, V = symbols('x y z V', real=True)
    
    def s3p(v1, v2, v3):
        "Scalar triple product of three vectors = v1⋅(v2 x v3)"
        return (v1[0] * (v2[1]*v3[2] - v2[2]*v3[1]) +
                v1[1] * (v2[2]*v3[0] - v2[0]*v3[2]) +
                v1[2] * (v2[0]*v3[1] - v2[1]*v3[0]))
    
    # Define vectors between vertices
    v1p = [x - x1, y - y1, z - z1]
    v2p = [x - x2, y - y2, z - z2]
    v12 = [x2 - x1, y2 - y1, z2 - z1]
    v13 = [x3 - x1, y3 - y1, z3 - z1]
    v14 = [x4 - x1, y4 - y1, z4 - z1]
    v23 = [x3 - x2, y3 - y2, z3 - z2]
    v24 = [x4 - x2, y4 - y2, z4 - z2]
    
    # Volumes, see http://www.cdsimpson.net/2014/10/barycentric-coordinates.html
    V1 = abs(s3p(v2p, v24, v23))/S(6)
    V2 = abs(s3p(v1p, v13, v14))/S(6)
    V3 = abs(s3p(v1p, v14, v12))/S(6)
    V4 = abs(s3p(v1p, v12, v13))/S(6)
    Vf = abs(s3p(v12, v13, v14))/S(6)
    
    # Barycentric coordinates from volumes
    L1f, L2f, L3f, L4f = V1/V, V2/V, V3/V, V4/V
    #L4f = 1 - L1f - L2f - L3f 
    L1, L2, L3, L4 = Symbol('L1')(x, y, z), Symbol('L2')(x, y, z), Symbol('L3')(x, y, z), Symbol('L4')(x, y, z)
    
    replacements = [(L1.diff(x),  L1f.diff(x)),
                    (L2.diff(x),  L2f.diff(x)),
                    (L3.diff(x),  L3f.diff(x)),
                    (L4.diff(x),  L4f.diff(x)),
                    (L1.diff(y),  L1f.diff(y)),
                    (L2.diff(y),  L2f.diff(y)),
                    (L3.diff(y),  L3f.diff(y)),
                    (L4.diff(y),  L4f.diff(y)),
                    (L1.diff(z),  L1f.diff(z)),
                    (L2.diff(z),  L2f.diff(z)),
                    (L3.diff(z),  L3f.diff(z)),
                    (L4.diff(z),  L4f.diff(z)),
                    (L1, 1/S(4)),
                    (L2, 1/S(4)),
                    (L3, 1/S(4)),
                    (L4, 1/S(4)),
                    (x, (x1 + x2 + x3 + x4)/4),
                    (y, (y1 + y2 + y3 + y4)/4),
                    (z, (z1 + z2 + z3 + z4)/4)]
    
    if degree == 1:
        # Barycentric linear shape functions on a triangle
        phi = [L1, L2, L3, L4]
        
        print('V = %s' % Vf)
        
        def print_code(name, func, index):
            code = ['# %s' % name]
            for i in range(4):
                v = func(phi[i])
                v = v.subs(replacements)
                v = v.simplify()
                code.append('A[icell,%d,%d] = %s' % (index, i, v))
            code.append('')
            print('\n'.join(code))
        
        print_code('Value at xc (also the cell average value)', lambda q: q, 0)
        print_code('d/dx', lambda q: q.diff(x), 1)
        print_code('d/dy', lambda q: q.diff(y), 2)
        print_code('d/dz', lambda q: q.diff(z), 3)
    
    elif degree == 2:
        


def produce_code_with_sympy_DG2_2D():
    """
    Use sympy to derive the coefficients needed to convert from DG2 polynomials
    on a triangle to a Taylor basis by expressing the derivatives at the centre
    point of the triangle
    (Only used for generating the code above, not used while running Ocellaris)
    """
    from sympy import Symbol, symbols, S
    x1, x2, x3, y1, y2, y3 = symbols('x1 x2 x3 y1 y2 y3')
    x, y, D = symbols('x y D')
    
    # Barycentric coordinates from cartesian corner coordinates
    L1f = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))/D
    L2f = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))/D
    L3f = 1 - L1f - L2f
    L1, L2, L3 = Symbol('L1')(x, y), Symbol('L2')(x, y), Symbol('L3')(x, y)
    
    replacements = [(L1.diff(x, x),  0),
                    (L2.diff(x, x),  0),
                    (L3.diff(x, x),  0),
                    (L1.diff(y, y),  0),
                    (L2.diff(y, y),  0),
                    (L3.diff(y, y),  0),
                    (L1.diff(x, y),  0),
                    (L2.diff(x, y),  0),
                    (L3.diff(x, y),  0),
                    (L1.diff(x),  L1f.diff(x)),
                    (L2.diff(x),  L2f.diff(x)),
                    (L3.diff(x),  L3f.diff(x)),
                    (L1.diff(y),  L1f.diff(y)),
                    (L2.diff(y),  L2f.diff(y)),
                    (L3.diff(y),  L3f.diff(y)),
                    (L1, 1/S(3)),
                    (L2, 1/S(3)),
                    (L3, 1/S(3))]
    
    # Barycentric quadratic shape functions on a triangle
    phi = [None]*6
    phi[0] = 2*L1*(L1 - 1/S(2))
    phi[1] = 2*L2*(L2 - 1/S(2))
    phi[2] = 2*L3*(L3 - 1/S(2))
    phi[3] = 4*L2*L3
    phi[4] = 4*L3*L1
    phi[5] = 4*L1*L2
    
    def print_code(name, func, index):
        code = ['# %s' % name]
        for i in range(6):
            v = func(phi[i])
            v = v.subs(replacements)
            v = v.simplify()
            code.append('A[icell,%d,%d] = %s' % (index, i, v))
        code.append('')
        print('\n'.join(code))
    
    print_code('Value at xc (MODIFY ME TO GET THE CELL AVERAGE)', lambda q: q, 0)
    print_code('d/dx', lambda q: q.diff(x), 1)
    print_code('d/dy', lambda q: q.diff(y), 2)
    print_code('d/dx^2', lambda q: q.diff(x, x), 3)
    print_code('d/dy^2', lambda q: q.diff(y, y), 4)
    print_code('d/dx*dy', lambda q: q.diff(x, y), 5)
    
    # Cell average terms
    BT_QUADRATURE_TABLE = {}
    BT_QUADRATURE_TABLE[1] = ([[1/S(3), 1/S(3)]],
                              [1.0])
    BT_QUADRATURE_TABLE[2] = ([[2/S(3), 1/S(6)], [1/S(6), 2/S(3)], [1/S(6), 1/S(6)]],
                              [1/S(3), 1/S(3), 1/S(3)])
    BT_QUADRATURE_TABLE[3] = ([[1/S(3), 1/S(3)], [S(3)/5, S(1)/5], [S(1)/5, S(3)/5], [S(1)/5, S(1)/5]],
                              [-S(9)/16, S(25)/48, S(25)/48, S(25)/48])
    
    xv = x1*L1 + x2*L2 + x3*L3
    yv = y1*L1 + y2*L2 + y3*L3
    xc = xv.subs({L1: 1/S(3), L2: 1/S(3), L3: 1/S(3)})
    yc = yv.subs({L1: 1/S(3), L2: 1/S(3), L3: 1/S(3)})
    
    exprs = [('bar_xx', (xv-xc)**2/2),
             ('bar_yy', (yv-yc)**2/2),
             ('bar_xy', (xv-xc)*(yv-yc))]
    for i in range(6):
        exprs.append(('A[icell,0,%d]' % i, phi[i]))
    
    order = 3
    points, weights = BT_QUADRATURE_TABLE[order]
    for name, expr in exprs:
        v = 0
        for point, weight in zip(points, weights):
            loc = {L1: point[0],
                   L2: point[1],
                   L3: 1 - point[0] - point[1]}
            v += weight*expr.subs(loc)
        v = v.simplify().simplify()
        v = v.subs({x1: Symbol('x1'), x2: Symbol('x2'), x3: Symbol('x3'),
                    y1: Symbol('y1'), y2: Symbol('y2'), y3: Symbol('y3')})
        print('%s = %s' % (name, v))


if __name__ == '__main__':
    funcs = (#produce_code_with_sympy_DG1_2D,
             #produce_code_with_sympy_DG2_2D,
             #produce_code_with_sympy_DG1_3D,
             produce_code_with_sympy_DG2_3D,)
    for func in funcs:
        print('='*80)
        print(func.__name__ + '\n' + '-'*40 + '\n')
        func()
        print()
