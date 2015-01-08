"""
The HRIC upwind/downwind blending sheme
"""
import numpy
import dolfin
from . import ConvectionScheme, register_convection_scheme

@register_convection_scheme('HRIC')
class ConvectionSchemeHric2D(ConvectionScheme):
    def update(self, t, dt, velocity):
        """
        Update the values of c at the faces, i.e the value
        that will be advected
        """
        a_cell_vec = self.alpha_function.vector()
        beta = self.blending_function
        beta_vec = beta.vector()
        
        # Reconstruct the gradient to calculate upstream values
        self.gradient_reconstructor.reconstruct()
        gradient = self.gradient_reconstructor.gradient
        neighbour_minval = self.gradient_reconstructor.neighbour_minval
        neighbour_maxval = self.gradient_reconstructor.neighbour_maxval
        
        Co_max = 0
        for facet in dolfin.facets(self.mesh):
            i = facet.index()
            
            # Find the local cells (the two cells sharing this face)
            connected_cells = self.con12(i)

            if len(connected_cells) != 2:
                # This should be an exterior facet (on ds)
                assert facet.exterior()
                beta_vec[self.dofmap[i]] = 0.0
                continue
            
            # Indices of the two local cells
            ic0, ic1 = connected_cells

            # Facet midpoint
            face_mp = self.facet_centroids[facet.index()]
            
            # Velocity at the midpoint (do not care which side of the face)
            ump = numpy.zeros(2, float)
            velocity.eval(ump, face_mp)

            # Midpoint of local cells
            cell0_mp = self.centroids[ic0]
            cell1_mp = self.centroids[ic1]
            mp_dist = cell1_mp - cell0_mp

            # Vector from local cell midpoint to face midpoint
            vec0 = face_mp - cell0_mp

            # Find a normal pointing out from local cell 0
            normalpt = facet.normal()
            normal = numpy.array([normalpt.x(), normalpt.y()], float)
            if numpy.dot(vec0, normal) < 0:
                normal *= -1

            # Find indices of downstream ("D") cell and central ("C") cell
            if numpy.dot(normal, ump) > 0:
                iaC = ic0
                iaD = ic1
                vec_to_upstream = -mp_dist
            else:
                iaC = ic1
                iaD = ic0
                vec_to_upstream = mp_dist
            
            # Find alpha in D and C cells
            aD = a_cell_vec[self.alpha_dofmap[iaD]]
            aC = a_cell_vec[self.alpha_dofmap[iaC]]
            
            # Extrapolate to the upwind U cell using the gradient
            # and ensure boundedness using the neighbouring values
            aU = aC + numpy.dot(gradient[iaC], vec_to_upstream)
            aU = max(aU, neighbour_minval[iaC])
            aU = min(aU, neighbour_maxval[iaC])
            
            # Calculate the Courant number
            dx = (mp_dist[0]**2 + mp_dist[1]**2)**0.5
            u = (ump[0]**2 + ump[1]**2)**0.5
            Co = u*dt/dx
            Co_max = max(Co_max, Co)
            
            if aD == aU or self.force_upwind:
                # No change in this area, use upstream value
                beta_vec[self.dofmap[i]] = 0.0
                continue
            
            # Introduce normalized variables
            tilde_aC = (aC - aU)/(aD - aU)
            
            if tilde_aC <= 0 or tilde_aC >= 1:
                # Only upwind is stable
                beta_vec[self.dofmap[i]] = 0.0
                continue
            elif 0 <= tilde_aC <= 0.5:
                # Blend upwind and downwind
                tilde_aF = 2*tilde_aC
            else:
                # Downwind
                tilde_aF = 1

            tilde_aF_star2 = tilde_aF

            # Calculate the downstream blending factor (0=upstream, 1=downstream)
            tilde_beta = (tilde_aF_star2 - tilde_aC)/(1 - tilde_aC)
            assert 0.0 <= tilde_beta <= 1.0
            beta_vec[self.dofmap[i]] = tilde_beta
        
        beta.vector()[:] = beta_vec
        print 'HRIC alpha_face  %10.5f %10.5f,  Co_max = %.3f' % (beta_vec.min(), beta_vec.max(), Co_max)
