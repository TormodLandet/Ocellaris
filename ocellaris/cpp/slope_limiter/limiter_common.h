#ifndef __SLOPE_LIMITER_COMMON_H
#define __SLOPE_LIMITER_COMMON_H

#include <cstdint>
#include <limits>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Core>


namespace dolfin
{


using IntVecIn = Eigen::Ref<const Eigen::VectorXi>;
using DoubleVecIn = Eigen::Ref<const Eigen::VectorXd>;


enum BoundaryDofType {
  NOT_ON_BOUNDARY = 0,
  DIRICHLET = 1,
  NEUMANN = 2,
  OTHER = 3
};


struct SlopeLimiterInput
{
  // --------------------------------------------------------------------------
  // Connectivity and dofs
  // --------------------------------------------------------------------------

  // Dimensions of the arrays
  int num_cells_owned;
  int max_neighbours;

  // Number of neighbours for each dof (dimension Ndof)
  std::vector<int> num_neighbours;

  // Neighbours for each dof (dimension Ndof * max_neighbours)
  std::vector<int> neighbours;

  // Dofs for each cell (dimension num_cells_owned * ndofs per cell)
  std::vector<int> cell_dofs;
  std::vector<int> cell_dofs_dg0;

  // Coordinates of the cell vertices
  std::vector<double> vertex_coords;

  // We can clamp the limited values to a given range
  double global_min = std::numeric_limits<double>::lowest();
  double global_max = std::numeric_limits<double>::max();

  void set_arrays(const int num_cells_owned,
                  const int max_neighbours,
                  IntVecIn num_neighbours,
                  IntVecIn neighbours,
                  IntVecIn cell_dofs,
                  IntVecIn cell_dofs_dg0,
                  DoubleVecIn vertex_coords)
  {
    this->num_cells_owned = num_cells_owned;
    this->max_neighbours = max_neighbours;

    this->num_neighbours.resize(num_neighbours.size());
    for (int i = 0; i < num_neighbours.size(); i++)
      this->num_neighbours[i] = num_neighbours[i];

    this->neighbours.resize(neighbours.size());
    for (int i = 0; i < neighbours.size(); i++)
      this->neighbours[i] = neighbours[i];

    this->cell_dofs.resize(cell_dofs.size());
    for (int i = 0; i < cell_dofs.size(); i++)
      this->cell_dofs[i] = cell_dofs[i];

    this->cell_dofs_dg0.resize(cell_dofs_dg0.size());
    for (int i = 0; i < cell_dofs_dg0.size(); i++)
      this->cell_dofs_dg0[i] = cell_dofs_dg0[i];

    this->vertex_coords.resize(vertex_coords.size());
    for (int i = 0; i < vertex_coords.size(); i++)
      this->vertex_coords[i] = vertex_coords[i];

    this->limit_cell.resize(limit_cell.size());
    for (unsigned int i = 0; i < limit_cell.size(); i++)
      this->limit_cell[i] = static_cast<std::int8_t>(limit_cell[i]);
  }

  // Should we limit a given cell. Look up with cell number and get 1 or 0
  std::vector<std::int8_t> limit_cell;

  void set_limit_cell(IntVecIn limit_cell)
  {
    this->limit_cell.resize(limit_cell.size());
    for (int i = 0; i < limit_cell.size(); i++)
      this->limit_cell[i] = static_cast<std::int8_t>(limit_cell[i]);
  }

  // --------------------------------------------------------------------------
  // Boundary conditions
  // --------------------------------------------------------------------------

  std::vector<BoundaryDofType> boundary_dof_type;
  std::vector<float> boundary_dof_value;
  bool enforce_boundary_conditions;

  void set_boundary_values(IntVecIn boundary_dof_type,
                           DoubleVecIn boundary_dof_value,
                           const bool enforce_bcs)
  {
    this->boundary_dof_type.resize(boundary_dof_type.size());
    for (int i = 0; i < boundary_dof_type.size(); i++)
      this->boundary_dof_type[i] = static_cast<BoundaryDofType>(boundary_dof_type[i]);

    this->boundary_dof_value.resize(boundary_dof_value.size());
    for (int i = 0; i < boundary_dof_value.size(); i++)
      this->boundary_dof_value[i] = boundary_dof_value[i];

    this->enforce_boundary_conditions = enforce_bcs;
  }
};



} // end namespace dolfin

#endif
