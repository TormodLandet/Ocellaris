#ifndef __GRADIENT_RECONSTRUCTION_H
#define __GRADIENT_RECONSTRUCTION_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/GenericDofMap.h>
//#include <fstream>
//#include <boost/lexical_cast.hpp>

namespace dolfin
{

void reconstruct_gradient(const Function& alpha_function,
                          const Array<int>& num_neighbours,
                          const int max_neighbours,
                          const Array<int>& neighbours,
                          const Array<double>& lstsq_matrices,
                          const Array<double>& lstsq_inv_matrices,
                          Function& gradient)
{
	const FunctionSpace& V = *alpha_function.function_space();
	const FunctionSpace& Vvec = *gradient.function_space();

	const Mesh& mesh = *V.mesh();
	const std::size_t ndim = mesh.geometry().dim();

	// Get dofmaps
	const std::vector<la_index> alpha_dofmap = V.dofmap()->dofs();
	std::vector<la_index> gradient_dofmap[ndim];
	for (int d = 0; d < ndim; d++)
	{
		gradient_dofmap[d] = Vvec[d]->dofmap()->dofs();
	}
	std::shared_ptr<const GenericVector> a_cell_vec = alpha_function.vector();
	std::shared_ptr<GenericVector> gradient_vec = gradient.vector();

	double ATdotB[ndim];
	double grad[ndim];
	int i = 0;
	for (CellIterator cell(mesh); !cell.end(); ++cell)
	{
		// Reset ATdotB
		for (int d = 0; d < ndim; d++)
		{
			ATdotB[d] = 0.0;
		}

		// Get the value in this cell
		const la_index idx = cell->index();
		const la_index dix = alpha_dofmap[idx];
		double a0 = (*a_cell_vec)[dix];

		// Compute the transpose(A)*B  matrix vector product
		const int Nnbs = num_neighbours[i];
		const int ncells = num_neighbours.size();
		int start = i*ndim*max_neighbours;
		for (int n = 0; n < Nnbs; n++)
		{
			const la_index nidx = neighbours[i*max_neighbours+n];
			const la_index ndix = alpha_dofmap[nidx];
			double aN = (*a_cell_vec)[ndix];
			for (int d = 0; d < ndim; d++)
			{
				ATdotB[d] += lstsq_matrices[start+d*max_neighbours+n]*(aN - a0);
			}
		}

		// Compute the inv(AT*A) * ATdotB matrix vector product
		start = i*ndim*ndim;
		for (int d = 0; d < ndim; d++)
		{
			grad[d] = 0.0;
			for (int d2 = 0; d2 < ndim; d2++)
			{
				grad[d] += lstsq_inv_matrices[start+d*ndim+d2]*ATdotB[d2];
			}
			const la_index didx2 = gradient_dofmap[d][idx];
			gradient_vec->set_local(&grad[d], 1, &didx2);
		}
		i++;
	}
	gradient_vec->apply("insert");
}

}

#endif
