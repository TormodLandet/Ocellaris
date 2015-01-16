#ifndef __GRADIENT_RECONSTRUCTION_H
#define __GRADIENT_RECONSTRUCTION_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
//#include <fstream>
//#include <boost/lexical_cast.hpp>

namespace dolfin
{
	void reconstruct_gradient(const Function& alpha_function,
							  const Array<int>& num_neighbours,
							  const Array<int>& neighbours,
							  const Array<double>& lstsq_matrices,
							  const Array<double>& lstsq_inv_matrices,
							  const Function& neighbour_minval,
							  const Function& neighbour_maxval,
							  const Function& gradient);
}

#endif
