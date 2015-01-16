#include "gradient_reconstruction.h"

namespace dolfin
{
void reconstruct_gradient(const Function& alpha_function,
                          const Array<int>& num_neighbours,
                          const Array<int>& neighbours,
                          const Array<double>& lstsq_matrices,
                          const Array<double>& lstsq_inv_matrices,
                          const Function& neighbour_minval,
                          const Function& neighbour_maxval,
                          const Function& gradient)
{
	const FunctionSpace& V = *alpha_function.function_space();
	const Mesh& mesh = *V.mesh();
	const std::size_t gdim = mesh.geometry().dim();
	std::shared_ptr<const GenericVector> a_cell_vec = alpha_function.vector();



 /*
  alpha_dofmap = V.dofmap().dofs()
  Vvec = gradient.function_space()
  gradient_dofmap0 = Vvec.sub(0).dofmap().dofs()
  gradient_dofmap1 = Vvec.sub(1).dofmap().dofs()

  np_gradient = gradient.vector().array()
  np_minvals = neighbour_minval.vector().array()
  np_maxvals = neighbour_maxval.vector().array()

  for i, cell in enumerate(dolfin.cells(mesh)):
      idx = cell.index()
      dix = alpha_dofmap[idx]
      Nnbs = num_neighbours[i]
      nbs = neighbours[i,:Nnbs]

      # Get the matrices
      AT = lstsq_matrices[i,:,:Nnbs]
      ATAI = lstsq_inv_matrices[i]
      a0  = a_cell_vec[alpha_dofmap[idx]]
      b = [(a_cell_vec[alpha_dofmap[ni]] - a0) for ni in nbs]
      b = numpy.array(b, float)

      # Store min and max values which can be used to enforce convective boundedness
      np_minvals[dix] = b.min()
      np_maxvals[dix] = b.max()

      # Calculate the and store the gradient
      g = numpy.dot(ATAI, numpy.dot(AT, b))
      np_gradient[gradient_dofmap0[idx]] = g[0]
      np_gradient[gradient_dofmap1[idx]] = g[1]

  gradient.vector()[:] = np_gradient
  neighbour_minval.vector()[:] = np_minvals
  neighbour_maxval.vector()[:] = np_maxvals
*/
}
}
