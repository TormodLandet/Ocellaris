#ifndef __SLOPE_LIMITER_BASIC_H
#define __SLOPE_LIMITER_BASIC_H

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/GenericDofMap.h>

namespace dolfin
{


void hierarchical_taylor_slope_limiter_dg1(const Array<int>& num_neighbours,
                                           const int num_cells_owned,
                                           const int max_neighbours,
                                           const Array<int>& neighbours,
                                           const Array<int>& cell_dofs,
                                           const Array<int>& cell_dofs_dg0,
                                           const Array<double>& vertex_coords,
                                           double* taylor_arr,
                                           double* alpha_arr)
{
  const int vstride = 8;
  const int dstride = 3;

  // Loop over all cells that are owned by this process
  for (int ic = 0; ic < num_cells_owned; ic++)
  {
    double alpha = 1.0;

    // The cell centre is stored as vertex 4
    double cx = vertex_coords[ic*vstride + 3*2 + 0];
    double cy = vertex_coords[ic*vstride + 3*2 + 1];

    // Get the Taylor values for this cell
    double center_phi = taylor_arr[cell_dofs[ic * dstride + 0]];
    double center_phix = taylor_arr[cell_dofs[ic * dstride + 1]];
    double center_phiy = taylor_arr[cell_dofs[ic * dstride + 2]];

    for (int ivert = 0; ivert < 3; ivert++)
    {
      // Calculate the value of phi at the vertex
      double dx = vertex_coords[ic*vstride + ivert*2 + 0] - cx;
      double dy = vertex_coords[ic*vstride + ivert*2 + 1] - cy;
      double vertex_value = center_phi + center_phix * dx + center_phiy * dy;

      // Find highest and lowest value in the connected neighbour cells
      int dof = cell_dofs[ic * dstride + ivert];
      double lo = center_phi;
      double hi = center_phi;
      for (int inb = 0; inb < num_neighbours[dof]; inb++)
      {
        int nb = neighbours[dof * max_neighbours + inb];
        double nb_val = taylor_arr[cell_dofs[nb * dstride + 0]];
        lo = std::min(lo, nb_val);
        hi = std::max(hi, nb_val);
      }

      // Compute the slope limiter coefficient alpha
      double a = 1.0;
      if (vertex_value > center_phi)
      {
        a = (hi - center_phi) / (vertex_value - center_phi);
      }
      else if (vertex_value < center_phi)
      {
        a = (lo - center_phi) / (vertex_value - center_phi);
      }
      alpha = std::min(alpha, a);
    }

    // Slope limit this cell
    alpha_arr[cell_dofs_dg0[ic]] = alpha;
    taylor_arr[cell_dofs[ic * dstride + 1]] *= alpha;
    taylor_arr[cell_dofs[ic * dstride + 2]] *= alpha;
  }
}


void hierarchical_taylor_slope_limiter_dg2(const Array<int>& num_neighbours,
                                           const int num_cells_owned,
                                           const int max_neighbours,
                                           const Array<int>& neighbours,
                                           const Array<int>& cell_dofs,
                                           const Array<int>& cell_dofs_dg0,
                                           const Array<double>& vertex_coords,
                                           double* taylor_arr,
                                           double* alpha1_arr,
                                           double* alpha2_arr)
{
  const int vstride = 8;
  const int dstride = 6;

  // Loop over all cells that are owned by this process
  for (int ic = 0; ic < num_cells_owned; ic++)
  {
    double alpha[3] = {1.0, 1.0, 1.0};

    // The cell centre is stored as vertex 4
    double cx = vertex_coords[ic*vstride + 3*2 + 0];
    double cy = vertex_coords[ic*vstride + 3*2 + 1];

    // Get the Taylor values for this cell
    double center_phi = taylor_arr[cell_dofs[ic * dstride + 0]];
    double center_phix = taylor_arr[cell_dofs[ic * dstride + 1]];
    double center_phiy = taylor_arr[cell_dofs[ic * dstride + 2]];
    double center_phixx = taylor_arr[cell_dofs[ic * dstride + 3]];
    double center_phiyy = taylor_arr[cell_dofs[ic * dstride + 4]];
    double center_phixy = taylor_arr[cell_dofs[ic * dstride + 5]];

    for (int itaylor = 0; itaylor < 3; itaylor++)
    {
      for (int ivert = 0; ivert < 3; ivert++)
      {
        // Calculate the value of phi or its gradient at the vertex
        double dx = vertex_coords[ic*vstride + ivert*2 + 0] - cx;
        double dy = vertex_coords[ic*vstride + ivert*2 + 1] - cy;
        double base_value, vertex_value;
        if (itaylor == 0)
        {
          // Function value at the vertex (linear reconstruction)
          vertex_value = center_phi + center_phix * dx + center_phiy * dy;
          base_value = center_phi;
        }
        else if (itaylor == 1)
        {
          // Derivative in x direction at the vertex (linear reconstruction)
          vertex_value = center_phix + center_phixx * dx + center_phixy * dy;
          base_value = center_phix;
        }
        else if (itaylor == 2)
        {
          // Derivative in y direction at the vertex (linear reconstruction)
          vertex_value = center_phiy + center_phiyy * dy + center_phixy * dx;
          base_value = center_phiy;
        }

        // Find highest and lowest value in the connected neighbour cells
        int dof = cell_dofs[ic * dstride + ivert];
        double lo = base_value;
        double hi = base_value;
        for (int inb = 0; inb < num_neighbours[dof]; inb++)
        {
          int nb = neighbours[dof * max_neighbours + inb];
          double nb_val = taylor_arr[cell_dofs[nb * dstride + itaylor]];
          lo = std::min(lo, nb_val);
          hi = std::max(hi, nb_val);
        }

        // Compute the slope limiter coefficient alpha
        double a = 1.0;
        if (vertex_value > base_value)
        {
          a = (hi - base_value) / (vertex_value - base_value);
        }
        else if (vertex_value < base_value)
        {
          a = (lo - base_value) / (vertex_value - base_value);
        }
        alpha[itaylor] = std::min(alpha[itaylor], a);
      }
    }

    // Compute alphas by the hierarchical method
    double alpha2 = std::min(alpha[1], alpha[2]);
    double alpha1 = std::max(alpha[0], alpha2);

    // Slope limit this cell
    alpha1_arr[cell_dofs_dg0[ic]] = alpha1;
    alpha2_arr[cell_dofs_dg0[ic]] = alpha2;
    taylor_arr[cell_dofs[ic * dstride + 1]] *= alpha1;
    taylor_arr[cell_dofs[ic * dstride + 2]] *= alpha1;
    taylor_arr[cell_dofs[ic * dstride + 3]] *= alpha2;
    taylor_arr[cell_dofs[ic * dstride + 4]] *= alpha2;
    taylor_arr[cell_dofs[ic * dstride + 5]] *= alpha2;
  }
}

}

#endif
