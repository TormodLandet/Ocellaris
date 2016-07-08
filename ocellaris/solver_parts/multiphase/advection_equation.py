# enconding: utf8
from __future__ import division
import dolfin
from dolfin import dot, div, grad, jump, dx, dS


class AdvectionEquation(object):
    def __init__(self, simulation, Vc, cp, cpp, u_conv, beta, time_coeffs, dirichlet_bcs):
        """
        This class assembles the advection equation for a scalar function c 
        """
        self.simulation = simulation
        self.Vc = Vc
        self.cp = cp
        self.cpp = cpp
        self.u_conv = u_conv
        self.beta = beta
        self.time_coeffs = time_coeffs
        self.dirichlet_bcs = dirichlet_bcs
        
        # Discontinuous or continuous elements
        Vc_family = Vc.ufl_element().family()
        self.colour_is_discontinuous = (Vc_family == 'Discontinuous Lagrange')
        
        # Create UFL forms
        self.define_advection_equation()
    
    def define_advection_equation(self):
        """
        Setup the advection equation for the colour function
        
        This implementation assembles the full LHS and RHS each time they are needed
        """
        sim = self.simulation
        mesh = sim.data['mesh']
        n = dolfin.FacetNormal(mesh)
        
        # Trial and test functions
        Vc = self.Vc
        c = dolfin.TrialFunction(Vc)
        d = dolfin.TestFunction(Vc)
        
        c1, c2, c3 = self.time_coeffs
        dt = sim.data['dt']
        u_conv = self.u_conv
        
        if not self.colour_is_discontinuous:
            # Continous Galerkin implementation of the advection equation
            # FIXME: add stabilization
            eq = (c1*c + c2*self.cp + c3*self.cpp)/dt*d*dx + div(c*u_conv)*d*dx
        
        elif self.beta is not None:
            # Upstream and downstream normal velocities
            w_nU = (dot(u_conv, n) + abs(dot(u_conv, n)))/2
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2
            
            if self.beta is not None:
                # Define the blended flux
                # The blending factor beta is not DG, so beta('+') == beta('-')
                b = self.beta('+')
                flux = (1-b)*jump(c*w_nU) + b*jump(c*w_nD)
            else:
                flux = jump(c*w_nU)
            
            # Discontinuous Galerkin implementation of the advection equation 
            eq = (c1*c + c2*self.cp + c3*self.cpp)/dt*d*dx \
                 - dot(c*u_conv, grad(d))*dx \
                 + flux*jump(d)*dS
            
            # Enforce Dirichlet BCs weakly
            for dbc in self.dirichlet_bcs:
                eq += w_nD*dbc.func()*d*dbc.ds()
                eq += w_nU*c*d*dbc.ds()
        
        else:
            # Downstream normal velocities
            w_nD = (dot(u_conv, n) - abs(dot(u_conv, n)))/2
            
            # Discontinuous Galerkin implementation of the advection equation
            eq = (c1*c + c2*self.cp + c3*self.cpp)/dt*d*dx
            
            # Convection integrated by parts two times to bring back the original
            # div form (this means we must subtract and add all fluxes)
            eq += div(c*u_conv)*d*dx
            
            # Replace downwind flux with upwind flux on downwind internal facets
            eq -= jump(w_nD*d)*jump(c)*dS
            
            # Replace downwind flux with upwind BC flux on downwind external facets
            for dbc in self.dirichlet_bcs:
                # Subtract the "normal" downwind flux
                eq -= w_nD*c*d*dbc.ds()
                # Add the boundary value upwind flux
                eq += w_nD*dbc.func()*d*dbc.ds()
        
        a, L = dolfin.system(eq)
        self.form_lhs = a
        self.form_rhs = L
        self.tensor_lhs = None
        self.tensor_rhs = None
    
    def assemble_lhs(self):
        if self.tensor_lhs is None:
            self.tensor_lhs = dolfin.assemble(self.form_lhs)
        else:
            dolfin.assemble(self.form_lhs, tensor=self.tensor_lhs)
        return self.tensor_lhs

    def assemble_rhs(self):
        if self.tensor_rhs is None:
            self.tensor_rhs = dolfin.assemble(self.form_rhs)
        else:
            dolfin.assemble(self.form_rhs, tensor=self.tensor_rhs)
        return self.tensor_rhs
