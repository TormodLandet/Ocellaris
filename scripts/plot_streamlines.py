"""
Read an Ocellaris restart file and plot streamlines
"""
import sys, os
import dolfin as df
from matplotlib import pyplot


class StreamFunction(object):
    def __init__(self, u, boundary_is_streamline=False, degree=1):
        """
        Heavily based on
        https://github.com/mikaem/fenicstools/blob/master/fenicstools/Streamfunctions.py
        
        Stream function for a given general 2D velocity field.
        The boundary conditions are weakly imposed through the term
        
            inner(q, grad(psi)*n)*ds, 
        
        where grad(psi) = [-v, u] is set on all boundaries. 
        This should work for any collection of boundaries: 
        walls, inlets, outlets etc.    
        """
        Vu = u[0].function_space()
        mesh = Vu.mesh()
        
        # Check dimension
        if not mesh.geometry().dim() == 2:
            df.error("Stream-function can only be computed in 2D.")
    
        # Define the weak form 
        V = df.FunctionSpace(mesh, 'CG', degree)
        q = df.TestFunction(V)
        psi = df.TrialFunction(V)
        n = df.FacetNormal(mesh)
        a = df.dot(df.grad(q), df.grad(psi))*df.dx
        L = df.dot(q, df.curl(u))*df.dx 
        
        if boundary_is_streamline: 
            # Strongly set psi = 0 on entire domain boundary
            self.bcs = [df.DirichletBC(V, df.Constant(0), df.DomainBoundary())]
            self.normalize = False
        else:
            self.bcs = []
            self.normalize = True
            L = L + q*(n[1]*u[0] - n[0]*u[1])*df.ds
            
        # Create preconditioned iterative solver
        solver = df.PETScKrylovSolver('gmres', 'hypre_amg')
        solver.parameters['nonzero_initial_guess'] = True
        solver.parameters['relative_tolerance'] = 1e-10
        solver.parameters['absolute_tolerance'] = 1e-10
        
        # Store for later computation
        self.psi = df.Function(V)
        self.A = df.assemble(a)
        self.L = L
        self.mesh = mesh
        self.solver = solver
        self._triangulation = None
    
    def compute(self):
        """
        Compute the stream function
        """
        b = df.assemble(self.L)
        
        if self.normalize:
            df.normalize(b)
        
        for bc in self.bcs:
            bc.apply(self.A, b)
        self.solver.solve(self.A, self.psi.vector(), b)
        
        if self.normalize: 
            df.normalize(self.psi.vector())
    
        return self.psi
    
    def plot(self, mpl_ax, levels=50, lw=0.3, mesh_alpha=0, mesh_lw=0.2):
        """
        Plot the function on a matplotlib axes. Call .compute() first
        to calculate the stream function
        """
        if self._triangulation is None:
            from matplotlib.tri import Triangulation
            coords = self.mesh.coordinates()
            triangles = []
            for cell in df.cells(self.mesh):
                cell_vertices = cell.entities(0)
                triangles.append(cell_vertices)
            self._triangulation = Triangulation(coords[:,0], coords[:,1], triangles)

        if mesh_alpha > 0:
            mpl_ax.triplot(self._triangulation, color='#000000', alpha=mesh_alpha, lw=mesh_lw)
        
        Z = self.psi.compute_vertex_values()
        if all(Z == 0):
            return
        
        mpl_ax.tricontour(self._triangulation, Z, levels, colors='#0000AA',
                          linewidths=lw, linestyles='solid')


def load_simulation(h5_file_name):
    assert os.path.isfile(h5_file_name)
    h5 = df.HDF5File(df.MPI.comm_world, h5_file_name, 'r')
    assert h5.has_dataset('ocellaris')
    
    mesh = df.Mesh()
    h5.read(mesh, '/mesh', False)
    
    u0_signature = h5.attributes('/u0')['signature']
    print u0_signature 
    eu = eval('df.' + u0_signature.replace('triangle', 'df.triangle'))
    Vu = df.FunctionSpace(mesh, eu)
    
    u0 = df.Function(Vu)
    u1 = df.Function(Vu)
    
    h5.read(u0, '/u0')
    h5.read(u1, '/u1')
    
    inp = h5.attributes('/ocellaris')['input_file']
    time = h5.attributes('/ocellaris')['time']
    h5.close()
    
    res = {'u0': u0,
           'u1': u1,
           'input': inp,
           'file_name': h5_file_name,
           'time': time}
    
    return res


def plot_stream_lines(res):
    u = df.as_vector([res['u0'], res['u1']])
    V = u[0].function_space()
    
    sf = StreamFunction(u, degree=V.ufl_element().degree())
    sf.compute()
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    sf.plot(ax, levels=50, lw=0.5, mesh_alpha=0.5, mesh_lw=0.2)
    pyplot.show()


if __name__ == '__main__':
    h5_file_name = sys.argv[1]
    res = load_simulation(h5_file_name)
    plot_stream_lines(res)
