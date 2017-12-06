#ifndef __GRADIENT_RECONSTRUCTION_H
#define __GRADIENT_RECONSTRUCTION_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/GenericDofMap.h>
#include <Eigen/Core>

namespace dolfin
{

using IntVecIn = Eigen::Ref<const Eigen::VectorXi>;
using IntMatIn = Eigen::Ref<const Eigen::MatrixXi>;
using DoubleVecIn = Eigen::Ref<const Eigen::VectorXd>;

void reconstruct_gradient(const Function& alpha_function,
                          IntVecIn num_neighbours,
                          IntMatIn neighbours,
                          DoubleVecIn lstsq_matrices,
                          DoubleVecIn lstsq_inv_matrices,
                          std::vector<std::shared_ptr<Function>>& gradient)
{
  const FunctionSpace& V = *alpha_function.function_space();
  const FunctionSpace& Vgrad = *gradient[0]->function_space();

  const Mesh& mesh = *V.mesh();
  const std::size_t ndim = mesh.geometry().dim();

  const std::size_t num_cells_owned = neighbours.rows();
  const std::size_t num_neighbours_max = neighbours.cols();

  // Get dofmaps
  const std::vector<la_index> alpha_dofmap = V.dofmap()->dofs();
  const std::vector<la_index> gradient_dofmap = Vgrad.dofmap()->dofs();

  // Get vectors
  std::vector<double> a_vec;
  alpha_function.vector()->get_local(a_vec);
  std::vector<std::vector<double>> grad_vec(ndim);
  for (std::size_t d = 0; d < ndim; d++)
    gradient[d]->vector()->get_local(grad_vec[d]);

  double ATdotB[ndim];
  double grad[ndim];
  for (int icell = 0; icell < num_cells_owned; icell++)
  {
    // Reset ATdotB
    for (std::size_t d = 0; d < ndim; d++)
    {
      ATdotB[d] = 0.0;
    }

    // Get the value in this cell
    const la_index dix = alpha_dofmap[icell];
    const double a0 = a_vec[dix];

    // Compute the transpose(A)*B  matrix vector product
    int start = icell * ndim * num_neighbours_max;
    for (int n = 0; n < num_neighbours[icell]; n++)
    {
      const la_index nidx = neighbours(icell, n);
      const la_index ndix = alpha_dofmap[nidx];
      double aN = a_vec[ndix];
      for (std::size_t d = 0; d < ndim; d++)
      {
        ATdotB[d] += lstsq_matrices[start + d * num_neighbours_max + n] * (aN - a0);
      }
    }

    // Compute the inv(AT*A) * ATdotB matrix vector product
    start = icell * ndim * ndim;
    for (std::size_t d = 0; d < ndim; d++)
    {
      grad[d] = 0.0;
      for (std::size_t d2 = 0; d2 < ndim; d2++)
      {
        grad[d] += lstsq_inv_matrices[start + d * ndim + d2] * ATdotB[d2];
      }
      const la_index dof = gradient_dofmap[icell];
      grad_vec[d][dof] = grad[d];
    }
  }
  for (std::size_t d = 0; d < ndim; d++)
  {
    gradient[d]->vector()->set_local(grad_vec[d]);
    gradient[d]->vector()->apply("insert");
  }
}

}

#endif
