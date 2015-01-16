"""
The HRIC upwind/downwind blending sheme
"""
import numpy
import dolfin
from . import ConvectionScheme, register_convection_scheme

@register_convection_scheme('HRIC')
class ConvectionSchemeHric2D(ConvectionScheme):
    
    description = 'High Resolution Interface Capturing (HRIC, Modified HRIC and Refined HRIC)'
    
    def __init__(self, simulation, func_name):
        """
        Implementation of the HRIC VOF convection scheme
        
        HRIC:
          "A Two-Fluid Navie-Stokes Solver to Simulate Water Entry"
          Twenty-Second Symposium on Naval Hydrodynamics, 1999
          S. Muzaferija, M. Peric, P. Sames, t. Schellin
        
        MHRIC:
          "Fluent Theory Guide v15.0"
          2013
          ANSYS 
        
        RHRIC:
          "Free surface flow analysis based on improved HRIC VOF method"
          Journal of the Society of Naval Architects of Korea
          Vol. 47, No 3, pp. 279-290, June 2010
          Il-Ryong Park, Kwang-Soo Kim, Jin Kim and Suak-Ho Van
        """
        super(ConvectionSchemeHric2D, self).__init__(simulation, func_name)
        self.variant = simulation.input.get('convection', {}).get(func_name, {}).get('HRIC_version', 'HRIC')
    
    def update(self, t, dt, velocity):
        """
        Update the values of the blending function beta at the facets
        according to the HRIC algorithm. Several versions of HRIC
        are implemented
        """
        a_cell_vec = self.alpha_function.vector()
        beta = self.blending_function
        beta_vec = beta.vector()
        
        conFC = self.simulation.data['connectivity_FC']
        facet_info = self.simulation.data['facet_info']
        cell_info = self.simulation.data['cell_info']
        
        # Reconstruct the gradient to calculate upstream values
        #self.gradient_reconstructor.reconstruct()
        gradient = self.gradient_reconstructor.gradient
        gradient_dofmap0 = self.gradient_reconstructor.gradient_dofmap0
        gradient_dofmap1 = self.gradient_reconstructor.gradient_dofmap1
        neighbour_minval = self.gradient_reconstructor.neighbour_minval
        neighbour_maxval = self.gradient_reconstructor.neighbour_maxval
        
        EPS = 1e-6
        Co_max = 0
        for facet in dolfin.facets(self.mesh):
            fidx = facet.index()
            finfo = facet_info[fidx]
            
            # Find the local cells (the two cells sharing this face)
            connected_cells = conFC(fidx)
            
            if len(connected_cells) != 2:
                # This should be an exterior facet (on ds)
                assert facet.exterior()
                beta_vec[self.dofmap[fidx]] = 0.0
                continue
            
            # Indices of the two local cells
            ic0, ic1 = connected_cells
            
            # Hack to speed up calculations
            nmin0 = neighbour_minval.vector()[self.alpha_dofmap[ic0]]
            nmax0 = neighbour_maxval.vector()[self.alpha_dofmap[ic0]]
            nmin1 = neighbour_minval.vector()[self.alpha_dofmap[ic1]]
            nmax1 = neighbour_maxval.vector()[self.alpha_dofmap[ic1]]
            if (nmax0 - nmin0) < EPS and (nmax1 - nmin1) < EPS:
                beta_vec[self.dofmap[fidx]] = 0.0
                continue
            
            # Velocity at the midpoint (do not care which side of the face)
            ump = numpy.zeros(2, float)
            velocity.eval(ump, finfo.midpoint)

            # Midpoint of local cells
            cell0_mp = cell_info[ic0].midpoint
            cell1_mp = cell_info[ic1].midpoint
            mp_dist = cell1_mp - cell0_mp

            # Normal pointing out of cell 0
            normal = finfo.normal

            # Find indices of downstream ("D") cell and central ("C") cell
            uf = numpy.dot(normal, ump) 
            if uf > 0:
                iaC = ic0
                iaD = ic1
                vec_to_downstream = mp_dist
                #nminC, nmaxC = nmin0, nmax0
            else:
                iaC = ic1
                iaD = ic0
                vec_to_downstream = -mp_dist
                #nminC, nmaxC = nmin1, nmax1
            
            # Find alpha in D and C cells
            aD = a_cell_vec[self.alpha_dofmap[iaD]]
            aC = a_cell_vec[self.alpha_dofmap[iaC]]
            
            # Gradient
            gdofs  = (gradient_dofmap0, gradient_dofmap1)
            func2vec = lambda fun, i: numpy.array([fun.vector()[dm[i]] for dm in gdofs], float)  
            gC = func2vec(gradient, iaC)
            len_gC2 = numpy.dot(gC, gC)
            
            # Upstream value
            # See Ubbink's PhD (1997) equations 4.21 and 4.22
            aU = aD - 2*numpy.dot(gC, vec_to_downstream)
            aU = min(max(aU, 0.0), 1.0)
            
            # Calculate the facet Courant number
            Co = abs(uf)*dt*finfo.area/cell_info[iaC].volume
            Co_max = max(Co_max, Co)
            
            if abs(aC - aD) < EPS or abs(aU - aD) < EPS:
                # No change in this area, use upstream value
                beta_vec[self.dofmap[fidx]] = 0.0
                continue
            
            # Angle between face normal and surface normal
            len_normal2 = numpy.dot(normal, normal)
            cos_theta = numpy.dot(normal, gC) / (len_normal2*len_gC2)**0.5
            
            # Introduce normalized variables
            tilde_aC = (aC - aU)/(aD - aU)
            
            if tilde_aC <= 0 or tilde_aC >= 1:
                # Only upwind is stable
                beta_vec[self.dofmap[fidx]] = 0.0
                continue
            
            if self.variant == 'HRIC':
                # Compressive scheme
                tilde_aF = 2*tilde_aC if 0 <= tilde_aC <= 0.5 else 1
                
                # Correct tilde_aF to avoid aligning with interfaces
                t = abs(cos_theta)**0.5
                tilde_aF_star = tilde_aF*t + tilde_aC*(1-t)
                
                # Correct tilde_af_star for high Courant numbers
                if Co < 0.4:
                    tilde_aF_final = tilde_aF_star
                elif Co < 0.75:
                    tilde_aF_final = tilde_aC + (tilde_aF_star - tilde_aC)*(0.75 - Co)/(0.75 - 0.4)
                else:
                    tilde_aF_final = tilde_aC
            
            elif self.variant == 'MHRIC':
                # Compressive scheme
                tilde_aF = 2*tilde_aC if 0 <= tilde_aC <= 0.5 else 1
                
                # Less compressive scheme
                tilde_aF_ultimate_quickest = min((6*tilde_aC + 3)/8, tilde_aF)
                
                # Correct tilde_aF to avoid aligning with interfaces
                t = abs(cos_theta)**0.5
                tilde_aF_final = tilde_aF*t + tilde_aF_ultimate_quickest*(1-t)
            
            elif self.variant == 'RHRIC':
                # Compressive scheme
                tilde_aF_hyperc = min(tilde_aC/Co, 1)
                
                # Less compressive scheme
                tilde_aF_hric = min(tilde_aC*Co + 2*tilde_aC*(1-Co), tilde_aF_hyperc)
                
                # Correct tilde_aF to avoid aligning with interfaces 
                t = cos_theta**4
                tilde_aF_final = tilde_aF_hyperc*t + tilde_aF_hric*(1-t)
            
            # Avoid tilde_aF being slightly lower that tilde_aC due to
            # floating point errors, it must be greater or equal 
            if tilde_aC - EPS < tilde_aF_final < tilde_aC:
                tilde_aF_final = tilde_aC
            
            # Calculate the downstream blending factor (0=upstream, 1=downstream)
            tilde_beta = (tilde_aF_final - tilde_aC)/(1 - tilde_aC)
            
            if not (0.0 <= tilde_beta <= 1.0):
                print 'ERROR, tilde_beta %r is out of range [0, 1]' % tilde_beta
                print ' face normal: %r' % normal
                print ' surface gradient: %r' % gC
                print ' cos(theta): %r' % cos_theta
                print ' sqrt(abs(cos(theta))) %r' % t
                print ' tilde_aF_final %r' % tilde_aF_final
                print ' tilde_aC %r' % tilde_aC
                print ' aU %r, aC %r, aD %r' % (aU, aC, aD)
            
            assert 0.0 <= tilde_beta <= 1.0
            beta_vec[self.dofmap[fidx]] = tilde_beta
        
        beta.vector()[:] = beta_vec
        self.simulation.reporting.report_timestep_value('Cof_max', Co_max)
        #print 'HRIC alpha_face  %10.5f %10.5f,  Co_max = %.3f' % (beta_vec.min(), beta_vec.max(), Co_max)
