import numpy
import dolfin as df
from dolfin import dot, ds, dS, dx
from ocellaris.utils import verify_key, get_dof_neighbours
from ocellaris.utils import lagrange_to_taylor, taylor_to_lagrange
from ocellaris.solver_parts.slope_limiter.hierarchical_taylor import HierarchicalTaylorSlopeLimiter
from . import register_velocity_slope_limiter, VelocitySlopeLimiterBase


@register_velocity_slope_limiter('LeastSquares')
class LeastSquaresSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    def __init__(self, simulation, vel, vel_name, use_cpp=True):
        """
        Use a Hierachical Taylor slope limiter on each of the velocity
        components and then least squares to reestablish the facet
        fluxes to ensure continuity 
        """
        # Verify input
        V = vel[0].function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        loc = 'SlopeLimiterVelocity'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (2,), loc)
        verify_key('function shape', vel.ufl_shape, [(2,)], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        
        # Store input
        self.simulation = simulation
        self.vel = vel
        self.vel_name = vel_name
        self.degree = degree
        self.mesh = mesh
        self.use_cpp = use_cpp
        
        # Create slope limiters for the velocity components
        self.temp_vels = [df.Function(V), df.Function(V)]
        sl0 = HierarchicalTaylorSlopeLimiter(vel_name, self.temp_vels[0], use_cpp=use_cpp, output_name=vel_name+'0')
        sl1 = HierarchicalTaylorSlopeLimiter(vel_name, self.temp_vels[1], use_cpp=use_cpp, output_name=vel_name+'1')
        self.slope_limiters = [sl0, sl1]
        
        # Fast access to cell dofs
        dm = V.dofmap()
        indices = list(range(self.mesh.num_cells()))
        self.cell_dofs = [dm.cell_dofs(i) for i in indices]
        
        self.additional_plot_funcs = sum((sl.additional_plot_funcs for sl in self.slope_limiters), [])
        
        # Define the over determined form of the projection
        self._define_form()
    
    def _define_form(self):
        V = self.vel[0].function_space()
        mesh = V.mesh()
        
        # The mixed function space of the projection test functions (non-square system)
        k = 2
        e1 = df.FiniteElement('DGT', mesh.ufl_cell(), 0)
        e2 = df.VectorElement('DG', mesh.ufl_cell(), 0)
        e3 = df.VectorElement('DG', mesh.ufl_cell(), k)
        em = df.MixedElement([e1, e2, e3])
        W = df.FunctionSpace(mesh, em)
        Vvec = df.VectorFunctionSpace(mesh, 'DG', 2)
        
        v1, v2, v3 = df.TestFunctions(W)
        u = df.TrialFunction(Vvec)
        n = df.FacetNormal(mesh)
        
        # Original and limited velocities
        w = self.vel
        m = df.as_vector(self.temp_vels)
        
        # Equation 1 - flux through the sides should not be changed
        a = L = 0
        for R in '+-':
            a += dot(u(R), n(R))*v1(R)*dS
            L += dot(w(R), n(R))*v1(R)*dS
        a += dot(u, n)*v1*ds
        L += dot(w, n)*v1*ds
        
        # Equation 2 - cell averages should not be changed
        a += dot(u, v2)*dx
        L += dot(w, v2)*dx
        
        # Equation 3 - velocity should be limited
        a += dot(u, v3)*dx
        L += dot(m, v3)*dx
        
        self.form = a, L
    
    def run(self):
        """
        Perform slope limiting of the velocity field
        """
        timer = df.Timer('Ocellaris LeastSquaresSlopeLimiterVelocity')
        
        # Perform initial slope limiting
        for d in range(2):
            self.temp_vels[d].assign(self.vel[d])
            self.slope_limiters[0].run()
            
        V = self.vel[0].function_space()
        mesh = V.mesh()
        
        u0_vals = self.vel[0].vector().get_local()
        u1_vals = self.vel[1].vector().get_local()
        
        a, L = self.form
        for cell in df.cells(mesh, 'regular'):
            dofs = self.cell_dofs[cell.index()]
            A = df.assemble_local(a, cell)
            b = df.assemble_local(L, cell)
            u = numpy.linalg.lstsq(A, b)[0]
            
            u0_vals[dofs] = u[:6]
            u1_vals[dofs] = u[6:]
        
        self.vel[0].vector().set_local(u0_vals)
        self.vel[1].vector().set_local(u1_vals)
        self.vel[0].vector().apply('insert')
        self.vel[1].vector().apply('insert')
        
        timer.stop()

@register_velocity_slope_limiter('HierarchicalTaylor')
class HierarchicalTaylorSlopeLimiterVelocity(VelocitySlopeLimiterBase):
    def __init__(self, simulation, vel, vel_name, use_cpp=True):
        """
        Limit the slope of the given vector field to obtain boundedness
        """
        # Verify input
        V = vel[0].function_space()
        mesh = V.mesh()
        family = V.ufl_element().family()
        degree = V.ufl_element().degree()
        loc = 'SlopeLimiterVelocity'
        verify_key('slope limited function', family, ['Discontinuous Lagrange'], loc)
        verify_key('slope limited degree', degree, (2,), loc)
        verify_key('function shape', vel.ufl_shape, [(2,)], loc)
        verify_key('topological dimension', mesh.topology().dim(), [2], loc)
        
        # Store input
        self.simulation = simulation
        self.vel = vel
        self.vel_name = vel_name
        self.degree = degree
        self.mesh = mesh
        self.use_cpp = use_cpp
        
        # Alpha factors are secondary outputs
        V0 = df.FunctionSpace(self.mesh, 'DG', 0)
        self.alpha_funcs = []
        for i in range(degree):
            func = df.Function(V0)
            name = 'SlopeLimiterVelocityAlpha%d_%s' % (i+1, vel_name)
            func.rename(name, name)
            self.alpha_funcs.append(func)
        self.additional_plot_funcs = self.alpha_funcs
        
        # Intermediate DG Taylor function space
        self.taylor = [df.Function(V), df.Function(V)]
        
        # Find the neighbour cells for each dof
        num_neighbours, neighbours = get_dof_neighbours(V)
        self.num_neighbours = num_neighbours
        self.neighbours = neighbours
        
        # Fast access to cell dofs
        dm, dm0 = V.dofmap(), V0.dofmap()
        indices = list(range(self.mesh.num_cells()))
        self.cell_dofs_V = [tuple(dm.cell_dofs(i)) for i in indices]
        self.cell_dofs_V0 = [int(dm0.cell_dofs(i)) for i in indices]
        
        # Find vertices for each cell
        mesh.init(2, 0)
        connectivity_CV = mesh.topology()(2, 0)
        vertices = []
        for ic in range(self.mesh.num_cells()):
            vnbs = tuple(connectivity_CV(ic))
            vertices.append(vnbs)
        self.vertices = vertices
        self.vertex_coordinates = mesh.coordinates()
    
    def run(self):
        """
        Perform slope limiting of DG Lagrange functions
        """
        timer = df.Timer('Ocellaris SlopeLimiterVelocity')
        self._run_dg2()
        timer.stop()
    
    def _run_dg2(self):
        """
        Perform slope limiting of a DG2 function
        """
        # Update the Taylor function space with the new DG values
        lagrange_to_taylor(self.vel[0], self.taylor[0])
        lagrange_to_taylor(self.vel[1], self.taylor[1])
        taylor_vals_u0 = self.taylor[0].vector().get_local()
        taylor_vals_u1 = self.taylor[1].vector().get_local()
        
        # Slope limiter coefficients
        alphas1 = self.alpha_funcs[0].vector().get_local()
        alphas2 = self.alpha_funcs[1].vector().get_local()
        
        V = self.vel[0].function_space()
        mesh = V.mesh()
        tdim = mesh.topology().dim()
        num_cells_owned = mesh.topology().ghost_offset(tdim)
        
        # Slope limit one cell at a time
        for icell in range(num_cells_owned):
            dofs = self.cell_dofs_V[icell]
            
            # vertex coordinates and the cell center
            cell_vertices = [self.vertex_coordinates[iv] for iv in self.vertices[icell]]
            center_pos_x = (cell_vertices[0][0] + cell_vertices[1][0] + cell_vertices[2][0]) / 3
            center_pos_y = (cell_vertices[0][1] + cell_vertices[1][1] + cell_vertices[2][1]) / 3
            
            # Find the minimum slope limiter coefficients alpha for each velocity component
            alpha_u0 = [1.0] * 3
            alpha_u1 = [1.0] * 3
            for taylor_vals, alpha in zip((taylor_vals_u0, taylor_vals_u1),
                                          (alpha_u0, alpha_u1)):
                
                center_values = [taylor_vals[dof] for dof in dofs]
                (center_phi, center_phix, center_phiy, center_phixx, 
                    center_phiyy, center_phixy) = center_values
                
                for taylor_dof in (0, 1, 2): 
                    for ivert in range(3):
                        dof = dofs[ivert]
                        dx = cell_vertices[ivert][0] - center_pos_x
                        dy = cell_vertices[ivert][1] - center_pos_y
                        
                        # Find vertex neighbours minimum and maximum values
                        base_value = center_values[taylor_dof]
                        minval = maxval = base_value
                        for nb in self.neighbours[dof]:
                            nb_center_val_dof = self.cell_dofs_V[nb][taylor_dof]
                            nb_val = taylor_vals[nb_center_val_dof]
                            minval = min(minval, nb_val)
                            maxval = max(maxval, nb_val)
                        
                        # Compute vertex value
                        if taylor_dof == 0:
                            # Function value at the vertex (linear reconstruction)
                            vertex_value = center_phi + center_phix * dx + center_phiy * dy
                        elif taylor_dof == 1:
                            # Derivative in x direction at the vertex  (linear reconstruction)
                            vertex_value = center_phix + center_phixx * dx + center_phixy * dy
                        else:
                            # Derivative in y direction at the vertex  (linear reconstruction)
                            vertex_value = center_phiy + center_phiyy * dy + center_phixy * dx
                        
                        # Compute the slope limiter coefficient alpha
                        if vertex_value > base_value:
                            a = (maxval - base_value) / (vertex_value - base_value)
                        elif vertex_value < base_value:
                            a = (minval - base_value) / (vertex_value - base_value)
                        alpha[taylor_dof] = min(alpha[taylor_dof], a)
            
            alpha2 = min(min(alpha_u0[1], alpha_u0[2]),
                         min(alpha_u1[1], alpha_u1[2]))
            alpha1 = min(max(alpha_u0[0], alpha2),
                         max(alpha_u1[0], alpha2))
            
            dof_dg0 = self.cell_dofs_V0[icell]    
            alphas1[dof_dg0] = alpha1
            alphas2[dof_dg0] = alpha2
            
            for taylor_vals in (taylor_vals_u0, taylor_vals_u1):
                taylor_vals[dofs[1]] *= alpha1
                taylor_vals[dofs[2]] *= alpha1
                taylor_vals[dofs[3]] *= alpha2
                taylor_vals[dofs[4]] *= alpha2
                taylor_vals[dofs[5]] *= alpha2
        
        # Update the DG Lagrange functions with the limited DG Taylor values
        for taylor_vals, taylor, vel in zip((taylor_vals_u0, taylor_vals_u1),
                                            self.taylor, self.vel):
            taylor.vector().set_local(taylor_vals)
            taylor.vector().apply('insert')
            taylor_to_lagrange(taylor, vel)
        
        # Update alpha functions
        for alpha, alpha_vals in zip(self.alpha_funcs, (alphas1, alphas2)):
            alpha.vector().set_local(alpha_vals)
            alpha.vector().apply('insert')
