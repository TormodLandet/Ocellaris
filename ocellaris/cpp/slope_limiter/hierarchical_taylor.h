#ifndef __SLOPE_LIMITER_HT_H
#define __SLOPE_LIMITER_HT_H

#include <cstdint>
#include <vector>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/fem/GenericDofMap.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Core>
#include "limiter_common.h" // remove_in_jit


namespace dolfin
{

using DoubleVec = Eigen::Ref<Eigen::VectorXd>;

template <int Ndim>  // Ndim is 2 for 2D triangles and 3 for 3D tetrahedra
void hierarchical_taylor_slope_limiter_dg1(const SlopeLimiterInput& input,
                                           DoubleVec taylor_arr,
                                           DoubleVec taylor_arr_old,
                                           DoubleVec alpha_arr)
{
  static_assert(Ndim == 2 or Ndim == 3, "Only 2D and 3D supported");
  const int num_cells_owned = input.num_cells_owned;
  const int max_neighbours = input.max_neighbours;
  const std::vector<int>& num_neighbours = input.num_neighbours;
  const std::vector<int>& neighbours = input.neighbours;
  const std::vector<int>& cell_dofs = input.cell_dofs;
  const std::vector<int>& cell_dofs_dg0 = input.cell_dofs_dg0;
  const std::vector<double>& vertex_coords = input.vertex_coords;
  const std::vector<std::int8_t>& limit_cell = input.limit_cell;
  const std::vector<BoundaryDofType>& boundary_dof_type = input.boundary_dof_type;
  const std::vector<float>& boundary_dof_value = input.boundary_dof_value;
  const double global_min = input.global_min;
  const double global_max = input.global_max;

  const int nvert = Ndim + 1;
  const int vstride = (nvert + 1) * Ndim; // 8 for 2D, 15 for 3D
  const int dstride = nvert;

  double cx, cy, cz=0.0, dx, dy, dz=0.0;
  double center_phi, center_phix, center_phiy, center_phiz=0.0;

  // Loop over all cells that are owned by this process
  for (int ic = 0; ic < num_cells_owned; ic++)
  {
    double alpha = 1.0;

    // The cell centre is stored as vertex 4 (2D) or 5 (3D)
    cx = vertex_coords[ic*vstride + nvert*Ndim + 0];
    cy = vertex_coords[ic*vstride + nvert*Ndim + 1];
    if (Ndim == 3)
      cz = vertex_coords[ic*vstride + nvert*Ndim + 2];

    // Get the Taylor values for this cell
    center_phi = taylor_arr[cell_dofs[ic * dstride + 0]];
    center_phix = taylor_arr[cell_dofs[ic * dstride + 1]];
    center_phiy = taylor_arr[cell_dofs[ic * dstride + 2]];
    if (Ndim == 3)
      center_phiz = taylor_arr[cell_dofs[ic * dstride + 3]];

    const bool skip_this_cell = (limit_cell[ic] == 0);
    if (!skip_this_cell)
    {
      for (int ivert = 0; ivert < nvert; ivert++)
      {
        // Calculate the value of phi at the vertex
        dx = vertex_coords[ic*vstride + ivert*Ndim + 0] - cx;
        dy = vertex_coords[ic*vstride + ivert*Ndim + 1] - cy;
        dz = 0.0;
        if (Ndim == 3)
          dz = vertex_coords[ic*vstride + ivert*Ndim + 2] - cz;
        double vertex_value = center_phi + center_phix * dx + center_phiy * dy + center_phiz * dz;

        // Find highest and lowest value in the connected neighbour cells
        int dof = cell_dofs[ic * dstride + ivert];
        double lo = center_phi;
        double hi = center_phi;
        for (int inb = 0; inb < num_neighbours[dof]; inb++)
        {
          int nb = neighbours[dof * max_neighbours + inb];
          int nb_dof = cell_dofs[nb * dstride + 0];
          double nb_val = taylor_arr[nb_dof];
          lo = std::min(lo, nb_val);
          hi = std::max(hi, nb_val);

          nb_val = taylor_arr_old[nb_dof];
          lo = std::min(lo, nb_val);
          hi = std::max(hi, nb_val);
        }

        // Modify local bounds to incorporate the boundary conditions
        if (boundary_dof_type[dof] == BoundaryDofType::DIRICHLET)
        {
          double bc_value = boundary_dof_value[dof];
          lo = std::min(lo, bc_value);
          hi = std::max(hi, bc_value);
        }

        // Modify local bounds to incorporate the global bounds
        lo = std::max(lo, global_min);
        hi = std::min(hi, global_max);
        center_phi = std::max(center_phi, global_min);
        center_phi = std::min(center_phi, global_max);

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
    }
    else
    {
      alpha = 1.0;
    }

    // Slope limit this cell
    alpha_arr[cell_dofs_dg0[ic]] = alpha;
    taylor_arr[cell_dofs[ic * dstride + 0]] = center_phi;
    taylor_arr[cell_dofs[ic * dstride + 1]] *= alpha;
    taylor_arr[cell_dofs[ic * dstride + 2]] *= alpha;
    if (Ndim == 3)
      taylor_arr[cell_dofs[ic * dstride + 3]] *= alpha;
  }
}


template <int Ndim>  // Ndim is 2 for 2D triangles and 3 for 3D tetrahedra
void hierarchical_taylor_slope_limiter_dg2(const SlopeLimiterInput& input,
                                           DoubleVec taylor_arr,
                                           DoubleVec taylor_arr_old,
                                           DoubleVec alpha1_arr,
                                           DoubleVec alpha2_arr)
{
  static_assert(Ndim == 2 or Ndim == 3, "Only 2D and 3D supported");
  const int num_cells_owned = input.num_cells_owned;
  const int max_neighbours = input.max_neighbours;
  const std::vector<int>& num_neighbours = input.num_neighbours;
  const std::vector<int>& neighbours = input.neighbours;
  const std::vector<int>& cell_dofs = input.cell_dofs;
  const std::vector<int>& cell_dofs_dg0 = input.cell_dofs_dg0;
  const std::vector<double>& vertex_coords = input.vertex_coords;
  const std::vector<std::int8_t>& limit_cell = input.limit_cell;
  const std::vector<BoundaryDofType>& boundary_dof_type = input.boundary_dof_type;
  const std::vector<float>& boundary_dof_value = input.boundary_dof_value;
  const double global_min = input.global_min;
  const double global_max = input.global_max;

  const int nvert = Ndim + 1;
  const int vstride = (nvert + 1) * Ndim; // 8 for 2D, 15 for 3D
  const int dstride = Ndim == 2 ? 6 : 10;

  double cx, cy, cz=0.0, dx, dy, dz=0.0;
  double center_phi, center_phix, center_phiy, center_phiz=0.0;
  double center_phixx, center_phiyy, center_phizz=0.0;
  double center_phixy, center_phixz=0.0, center_phiyz=0.0;

  // Loop over all cells that are owned by this process
  for (int ic = 0; ic < num_cells_owned; ic++)
  {
    double alpha[3] = {1.0, 1.0, 1.0};

    // The cell centre is stored as vertex 4 (2D) or 5 (3D)
    cx = vertex_coords[ic*vstride + nvert*Ndim + 0];
    cy = vertex_coords[ic*vstride + nvert*Ndim + 1];
    if (Ndim == 3)
      cz = vertex_coords[ic*vstride + nvert*Ndim + 2];

    // Get the Taylor values for this cell
    if (Ndim == 2)
    {
      center_phi = taylor_arr[cell_dofs[ic * dstride + 0]];
      center_phix = taylor_arr[cell_dofs[ic * dstride + 1]];
      center_phiy = taylor_arr[cell_dofs[ic * dstride + 2]];
      center_phixx = taylor_arr[cell_dofs[ic * dstride + 3]];
      center_phiyy = taylor_arr[cell_dofs[ic * dstride + 4]];
      center_phixy = taylor_arr[cell_dofs[ic * dstride + 5]];
    }
    else
    {
      center_phi = taylor_arr[cell_dofs[ic * dstride + 0]];
      center_phix = taylor_arr[cell_dofs[ic * dstride + 1]];
      center_phiy = taylor_arr[cell_dofs[ic * dstride + 2]];
      center_phiz = taylor_arr[cell_dofs[ic * dstride + 3]];
      center_phixx = taylor_arr[cell_dofs[ic * dstride + 4]];
      center_phiyy = taylor_arr[cell_dofs[ic * dstride + 5]];
      center_phizz = taylor_arr[cell_dofs[ic * dstride + 6]];
      center_phixy = taylor_arr[cell_dofs[ic * dstride + 7]];
      center_phixz = taylor_arr[cell_dofs[ic * dstride + 8]];
      center_phiyz = taylor_arr[cell_dofs[ic * dstride + 9]];
    }

    const bool skip_this_cell = (limit_cell[ic] == 0);
    for (int itaylor = 0; itaylor < nvert; itaylor++)
    {
      if (skip_this_cell)
        break;

      for (int ivert = 0; ivert < nvert; ivert++)
      {
        // Calculate the value of phi or its gradient at the vertex
        dx = vertex_coords[ic*vstride + ivert*Ndim + 0] - cx;
        dy = vertex_coords[ic*vstride + ivert*Ndim + 1] - cy;
        dz = 0.0;
        if (Ndim ==3)
          dz = vertex_coords[ic*vstride + ivert*Ndim + 2] - cz;

        double base_value, vertex_value;
        if (itaylor == 0)
        {
          // Function value at the vertex (linear reconstruction)
          vertex_value = center_phi + center_phix * dx + center_phiy * dy + center_phiz * dz;
          base_value = center_phi;
        }
        else if (itaylor == 1)
        {
          // Derivative in x direction at the vertex (linear reconstruction)
          vertex_value = center_phix + center_phixx * dx + center_phixy * dy + center_phixz * dz;
          base_value = center_phix;
        }
        else if (itaylor == 2)
        {
          // Derivative in y direction at the vertex (linear reconstruction)
          vertex_value = center_phiy + center_phiyy * dy + center_phixy * dx + center_phiyz * dz;
          base_value = center_phiy;
        }
        else if (itaylor == 3)
        {
          // Derivative in z direction at the vertex (linear reconstruction)
          vertex_value = center_phiz + center_phizz * dz + center_phixz * dx + center_phiyz * dy;
          base_value = center_phiz;
        }

        // Find highest and lowest value in the connected neighbour cells
        int dof = cell_dofs[ic * dstride + ivert];
        double lo = base_value;
        double hi = base_value;
        for (int inb = 0; inb < num_neighbours[dof]; inb++)
        {
          int nb = neighbours[dof * max_neighbours + inb];
          int nb_dof = cell_dofs[nb * dstride + itaylor];
          double nb_val = taylor_arr[nb_dof];
          lo = std::min(lo, nb_val);
          hi = std::max(hi, nb_val);

          nb_val = taylor_arr_old[nb_dof];
          lo = std::min(lo, nb_val);
          hi = std::max(hi, nb_val);
        }

        // Handle boundary conditions and global bounds
        bool dof_is_dirichlet = boundary_dof_type[dof] == BoundaryDofType::DIRICHLET;
        if (itaylor == 0)
        {
          // Modify local bounds to incorporate the boundary conditions
          if (dof_is_dirichlet)
          {
            // Value in the center of a mirrored cell on the other side of the boundary
            double bc_value = 2*boundary_dof_value[dof] - center_phi;
            lo = std::min(lo, bc_value);
            hi = std::max(hi, bc_value);
          }

          // Modify local bounds to incorporate the global bounds
          lo = std::max(lo, global_min);
          hi = std::min(hi, global_max);
          center_phi = std::max(center_phi, global_min);
          center_phi = std::min(center_phi, global_max);
        }
        else if (itaylor == 1 && dof_is_dirichlet)
        {
          // The derivative in the x-direction at the center of a mirrored cell
          double ddx = (boundary_dof_value[dof] - center_phi) / dx;
          double ddx2 = 4*ddx - 3*center_phix;
          lo = std::min(lo, ddx2);
          hi = std::max(hi, ddx2);
        }
        else if (itaylor == 2 && dof_is_dirichlet)
        {
          // The derivative in the y-direction at the center of a mirrored cell
          double ddy = (boundary_dof_value[dof] - center_phi) / dy;
          double ddy2 = 4*ddy - 3*center_phiy;
          lo = std::min(lo, ddy2);
          hi = std::max(hi, ddy2);
        }
        else if (itaylor == 3 && dof_is_dirichlet)
        {
          // The derivative in the z-direction at the center of a mirrored cell
          double ddz = (boundary_dof_value[dof] - center_phi) / dz;
          double ddz2 = 4*ddz - 3*center_phiz;
          lo = std::min(lo, ddz2);
          hi = std::max(hi, ddz2);
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

    double alpha1 = 1.0, alpha2 = 1.0;
    if (!skip_this_cell)
    {
      // Compute alphas by the hierarchical method
      alpha2 = std::min(alpha[1], alpha[2]);
      alpha1 = std::max(alpha[0], alpha2);
    }

    // Slope limit this cell
    alpha1_arr[cell_dofs_dg0[ic]] = alpha1;
    alpha2_arr[cell_dofs_dg0[ic]] = alpha2;

    if (Ndim == 2)
    {
      taylor_arr[cell_dofs[ic * dstride + 0]] = center_phi;
      taylor_arr[cell_dofs[ic * dstride + 1]] *= alpha1;
      taylor_arr[cell_dofs[ic * dstride + 2]] *= alpha1;
      taylor_arr[cell_dofs[ic * dstride + 3]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 4]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 5]] *= alpha2;
    }
    else
    {
      taylor_arr[cell_dofs[ic * dstride + 0]] = center_phi;
      taylor_arr[cell_dofs[ic * dstride + 1]] *= alpha1;
      taylor_arr[cell_dofs[ic * dstride + 2]] *= alpha1;
      taylor_arr[cell_dofs[ic * dstride + 3]] *= alpha1;
      taylor_arr[cell_dofs[ic * dstride + 4]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 5]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 6]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 7]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 8]] *= alpha2;
      taylor_arr[cell_dofs[ic * dstride + 9]] *= alpha2;
    }
  }
}


PYBIND11_MODULE(SIGNATURE, m)
{
  pybind11::class_<SlopeLimiterInput>(m, "SlopeLimiterInput")
      .def(pybind11::init())
      .def_readwrite("global_min", &SlopeLimiterInput::global_min)
      .def_readwrite("global_max", &SlopeLimiterInput::global_max)
      .def("set_arrays", &SlopeLimiterInput::set_arrays)
      .def("set_limit_cell", &SlopeLimiterInput::set_limit_cell)
      .def("set_boundary_values", &SlopeLimiterInput::set_boundary_values);
  m.def("hierarchical_taylor_slope_limiter_dg1_2D", &hierarchical_taylor_slope_limiter_dg1<2>);
  m.def("hierarchical_taylor_slope_limiter_dg1_3D", &hierarchical_taylor_slope_limiter_dg1<3>);
  m.def("hierarchical_taylor_slope_limiter_dg2_2D", &hierarchical_taylor_slope_limiter_dg2<2>);
  m.def("hierarchical_taylor_slope_limiter_dg2_3D", &hierarchical_taylor_slope_limiter_dg2<3>);
}


}

#endif
