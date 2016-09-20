"""
Convert to and from a Taylor basis
"""
from __future__ import division
import numpy
from dolfin import cells


def DG2_to_taylor(u, t):
    """
    Convert a function u (which is a DG2 function in 2D) to a DG Taylor basis
    which is stored in the DG2 function t. Most operations on the function t
    will be wrong since dolfin does not support Taylor basis DG elements. The
    t function should hence be read and modified by specially crafted code only! 
    """
    V = u.function_space()
    mesh = V.mesh()
    vertices = mesh.coordinates()
    dm = V.dofmap()
    dof_vals = u.vector().get_local()
    dof_vals_taylor = numpy.zeros_like(dof_vals)
    
    taylor_vals = numpy.zeros(6, float)
    
    for cell in cells(mesh):
        dofs = dm.cell_dofs(cell.index())
        vals = dof_vals[dofs]
        
        verts = cell.entities(0)
        x = numpy.zeros((6, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = (x[1] +  x[2])/2
        x[4] = (x[0] +  x[2])/2
        x[5] = (x[0] +  x[1])/2
                
        ###############################3
        # From sympy code gen
        
        ((x1, y1), (x2, y2), (x3, y3)) = x[:3]
        D = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
                
        # Value at xc
        taylor_vals[0] = (-1/9)*vals[0]
        taylor_vals[0] += (-1/9)*vals[1]
        taylor_vals[0] += (-1/9)*vals[2]
        taylor_vals[0] += (4/9)*vals[3]
        taylor_vals[0] += (4/9)*vals[4]
        taylor_vals[0] += (4/9)*vals[5]
        
        # d/dx
        taylor_vals[1] = ((y2 - y3)/(3*D))*vals[0]
        taylor_vals[1] += ((-y1 + y3)/(3*D))*vals[1]
        taylor_vals[1] += ((y1 - y2)/(3*D))*vals[2]
        taylor_vals[1] += (4*(-y2 + y3)/(3*D))*vals[3]
        taylor_vals[1] += (4*(y1 - y3)/(3*D))*vals[4]
        taylor_vals[1] += (4*(-y1 + y2)/(3*D))*vals[5]
        
        # d/dy
        taylor_vals[2] = ((-x2 + x3)/(3*D))*vals[0]
        taylor_vals[2] += ((x1 - x3)/(3*D))*vals[1]
        taylor_vals[2] += ((-x1 + x2)/(3*D))*vals[2]
        taylor_vals[2] += (4*(x2 - x3)/(3*D))*vals[3]
        taylor_vals[2] += (4*(-x1 + x3)/(3*D))*vals[4]
        taylor_vals[2] += (4*(x1 - x2)/(3*D))*vals[5]
        
        # d/dx^2
        taylor_vals[3] = (4*(y2 - y3)**2/D**2)*vals[0]
        taylor_vals[3] += (4*(y1 - y3)**2/D**2)*vals[1]
        taylor_vals[3] += (4*(-y1 + y2)**2/D**2)*vals[2]
        taylor_vals[3] += (8*(-y1 + y2)*(y1 - y3)/D**2)*vals[3]
        taylor_vals[3] += (8*(y1 - y2)*(y2 - y3)/D**2)*vals[4]
        taylor_vals[3] += (-8*(y1 - y3)*(y2 - y3)/D**2)*vals[5]
        
        # d/dy^2
        taylor_vals[4] = (4*(x2 - x3)**2/D**2)*vals[0]
        taylor_vals[4] += (4*(x1 - x3)**2/D**2)*vals[1]
        taylor_vals[4] += (4*(x1 - x2)**2/D**2)*vals[2]
        taylor_vals[4] += (-8*(x1 - x2)*(x1 - x3)/D**2)*vals[3]
        taylor_vals[4] += (8*(x1 - x2)*(x2 - x3)/D**2)*vals[4]
        taylor_vals[4] += (-8*(x1 - x3)*(x2 - x3)/D**2)*vals[5]
        
        # d/dx*dy
        taylor_vals[5] = (-4*(x2 - x3)*(y2 - y3)/D**2)*vals[0]
        taylor_vals[5] += (-4*(x1 - x3)*(y1 - y3)/D**2)*vals[1]
        taylor_vals[5] += (4*(x1 - x2)*(-y1 + y2)/D**2)*vals[2]
        taylor_vals[5] += (4*((x1 - x2)*(y1 - y3) + (x1 - x3)*(y1 - y2))/D**2)*vals[3]
        taylor_vals[5] += (-(4*(x1 - x2)*(y2 - y3) + 4*(x2 - x3)*(y1 - y2))/D**2)*vals[4]
        taylor_vals[5] += (4*((x1 - x3)*(y2 - y3) + (x2 - x3)*(y1 - y3))/D**2)*vals[5]
        
        dof_vals_taylor[dofs] = taylor_vals
    
    t.vector().set_local(dof_vals_taylor)


def taylor_to_DG2(t, u):
    """
    Take a function t as produced by DG2_to_taylor and convert it back
    to a standard Lagrange DG2 function, u.
    """
    V = u.function_space()
    mesh = V.mesh()
    vertices = mesh.coordinates()
    dm = V.dofmap()
    dof_vals_taylor = t.vector().get_local()
    dof_vals = numpy.zeros_like(dof_vals_taylor)
    
    vals = numpy.zeros(6, float)
    
    for cell in cells(mesh):
        dofs = dm.cell_dofs(cell.index())
        taylor_vals = dof_vals_taylor[dofs]
        
        # Corner and center coordinates
        verts = cell.entities(0)
        x = numpy.zeros((6, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = (x[1] +  x[2])/2
        x[4] = (x[0] +  x[2])/2
        x[5] = (x[0] +  x[1])/2
        xc = (x[0] + x[1] +  x[2])/3
        
        # Evaluate the Taylor basis at each node
        for i in range(6):
            dx, dy = x[i,0] - xc[0], x[i,1] - xc[1]
            vals[i] = taylor_vals[0] + dx*taylor_vals[1] + dy*taylor_vals[2] \
                      + dx**2/2*taylor_vals[3] + dy**2/2*taylor_vals[4] \
                      + dx*dy*taylor_vals[5]
        dof_vals[dofs] = vals
    
    u2.vector().set_local(dof_vals)


#########################################################################################

def DG2_to_taylor_numpy(u, t):
    """
    Convert DG2 function on triangles to a Taylor basis via linear solves on each element
    (only used for verification of the above direct code)
    """
    V = u.function_space()
    mesh = V.mesh()
    vertices = mesh.coordinates()
    dm = V.dofmap()
    dof_vals = u.vector().get_local()
    dof_vals_taylor = numpy.zeros_like(dof_vals)
    
    A = numpy.zeros((6, 6), float)
    
    for cell in cells(mesh):
        dofs = dm.cell_dofs(cell.index())
        vals = dof_vals[dofs]
        
        verts = cell.entities(0)
        x = numpy.zeros((6, 2), float)
        x[0] = vertices[verts[0]]
        x[1] = vertices[verts[1]]
        x[2] = vertices[verts[2]]
        x[3] = (x[1] +  x[2])/2
        x[4] = (x[0] +  x[2])/2
        x[5] = (x[0] +  x[1])/2
        xc = (x[0] + x[1] +  x[2])/3
        
        for i in range(6):
            dx, dy = x[i,0] - xc[0], x[i,1] - xc[1]
            A[i,0] = 1
            A[i,1] = dx
            A[i,2] = dy
            A[i,3] = dx**2/2
            A[i,4] = dy**2/2
            A[i,5] = dx*dy
        
        taylor_vals = numpy.linalg.solve(A, vals)
        dof_vals_taylor[dofs] = taylor_vals
    
    t.vector().set_local(dof_vals_taylor)


def produce_code_with_sympy():
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
            
            s = 'taylor_vals[%d]' % index
            s += ' = ' if i == 0 else ' += '
            s += '(%s)*vals[%d]' % (str(v), i)
            code.append(s)
        code.append('')
        print '\n'.join(code)
    
    print_code('Value at xc', lambda q: q, 0)
    print_code('d/dx', lambda q: q.diff(x), 1)
    print_code('d/dy', lambda q: q.diff(y), 2)
    print_code('d/dx^2', lambda q: q.diff(x, x), 3)
    print_code('d/dy^2', lambda q: q.diff(y, y), 4)
    print_code('d/dx*dy', lambda q: q.diff(x, y), 5)


if __name__ == '__main__':
    produce_code_with_sympy()
    print '-------------------------------------------------------------\n'
    
    from dolfin import UnitTriangleMesh, FunctionSpace, Function, errornorm
    
    mesh = UnitTriangleMesh()
    V = FunctionSpace(mesh, 'DG', 2)
    
    u, u2, t, tn = Function(V), Function(V), Function(V), Function(V)
    for vec in (numpy.array([0, 1, 2, 3, 4, 5], float),
                numpy.random.rand(6)): 
        u.vector().set_local(vec)
        
        # Convert to Taylor basis
        DG2_to_taylor_numpy(u, tn)
        DG2_to_taylor(u, t)
        print tn.vector().get_local(), 'numpy solve'
        print t.vector().get_local(), 'sympy expressions'
        print 'Error norm 1:', errornorm(tn, t, degree_rise=0)
        
        # Convert to back to DG basis
        taylor_to_DG2(t, u2)
        print u.vector().get_local(), 'original'
        print u2.vector().get_local(), 'round trip'
        print 'Error norm 2:', errornorm(u, u2, degree_rise=0)
        print
