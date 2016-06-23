#ifndef __SLOPE_LIMITER_BASIC_H
#define __SLOPE_LIMITER_BASIC_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/GenericDofMap.h>

namespace dolfin
{

void slope_limiter_basic_cg1(const Array<int>& num_neighbours,
                             const int num_cells_all,
                             const int num_cells_owned,
                             const int max_neighbours,
                             const Array<int>& neighbours,
                             const Array<int>& cell_dofs,
                             const Array<int>& cell_dofs_dg0,
                             const Array<int>& vertices,
                             double* exceedances,
                             double* results)
{
  // Calculate cell averages (which we will make sure to keep unchanged)
  std::vector<double> averages(num_cells_all);
  for (int ic = 0; ic < num_cells_all; ic++)
  {
    double avg = 0.0;
    for (int iv = 0; iv < 3; iv++)
    {
      avg += results[cell_dofs[ic*3 + iv]];
    }
    averages[ic] = avg/3.0;
  }

  // Modify vertex values
  for (int ic = 0; ic < num_cells_owned; ic++)
  {
    double avg = averages[ic];
    double exceedance = 0.0;
    double vals[3];

    // Check each vertex in this cell
    for (int iv = 0; iv < 3; iv++)
    {
      int vtx = vertices[ic*3 + iv];
      double vtx_val = results[cell_dofs[ic*3 + iv]];

      // Find highest and lowest value in the connected neighbour cells
      double lo = 1e100;
      double hi = -1e100;
      for (int inb = 0; inb < num_neighbours[vtx]; inb++)
      {
        int nb = neighbours[vtx*max_neighbours + inb];
        double nb_avg = averages[nb];
        lo = std::min(lo, nb_avg);
        hi = std::max(hi, nb_avg);
      }

      // Treat out of bounds values
      if (vtx_val < lo)
      {
        vals[iv] = lo;
        double ex = vtx_val - lo;
        if (std::abs(ex) > std::abs(exceedance)) exceedance = ex;
      }
      else if (vtx_val > hi)
      {
        vals[iv] = hi;
        double ex = vtx_val - hi;
        if (std::abs(ex) > std::abs(exceedance)) exceedance = ex;
      }
      else
      {
        vals[iv] = vtx_val;
      }
    }

    // Store the maximum absolute cell exceedance and quit early if possible
    exceedances[cell_dofs_dg0[ic]] = exceedance;
    if (exceedance == 0.0) continue;

    // Find the new average and which vertices can be adjusted to obtain the correct average
    double new_avg = (vals[0] + vals[1] + vals[2])/3.0;
    double eps = 0;
    bool moddable[3] = {false, false, false};
    if (std::abs(avg - new_avg) > 1e-15)
    {
     for(int iv = 0; iv < 3; iv++)
     {
       moddable[iv] = (new_avg > avg && vals[iv] > avg) || (new_avg < avg and vals[iv] < avg);
     }
     // Get number of vertex values that can be modified
     int nmod = ((int) moddable[0]) + ((int) moddable[1]) + ((int) moddable[2]);
     if (nmod == 0) {
       dolfin_assert(std::abs(exceedance) < 1e-14);
     }
     else
     {
       eps = (avg - new_avg)*3/nmod;
     }
    }

    // Update the result array
    for(int iv = 0; iv < 3; iv++)
    {
      results[cell_dofs[ic*3 + iv]] = vals[iv] + eps*((int)moddable[iv]);
    }
  }
}

}

#endif
