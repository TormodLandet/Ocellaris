#ifndef __SLOPE_LIMITER_COMMON_H
#define __SLOPE_LIMITER_COMMON_H

#include <limits>
#include <vector>


namespace dolfin
{


struct SlopeLimiterInput
{
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

  // Should we limit a given cell. Look up with cell number and get 1 or 0
  std::vector<short> limit_cell;

  // We can clamp the limited values to a given range
  double global_min = std::numeric_limits<double>::lowest();
  double global_max = std::numeric_limits<double>::max();

  void set_arrays(const int num_cells_owned,
                  const int max_neighbours,
                  const Array<int>& num_neighbours,
                  const Array<int>& neighbours,
                  const Array<int>& cell_dofs,
                  const Array<int>& cell_dofs_dg0,
                  const Array<double>& vertex_coords,
                  const Array<int>& limit_cell)
  {
    this->num_cells_owned = num_cells_owned;
    this->max_neighbours = max_neighbours;

    this->num_neighbours.reserve(num_neighbours.size());
    for (int i = 0; i < num_neighbours.size(); i++)
      this->num_neighbours.push_back(num_neighbours[i]);

    this->neighbours.reserve(neighbours.size());
    for (int i = 0; i < neighbours.size(); i++)
      this->neighbours.push_back(neighbours[i]);

    this->cell_dofs.reserve(cell_dofs.size());
    for (int i = 0; i < cell_dofs.size(); i++)
      this->cell_dofs.push_back(cell_dofs[i]);

    this->cell_dofs_dg0.reserve(cell_dofs_dg0.size());
    for (int i = 0; i < cell_dofs_dg0.size(); i++)
      this->cell_dofs_dg0.push_back(cell_dofs_dg0[i]);

    this->vertex_coords.reserve(vertex_coords.size());
    for (int i = 0; i < vertex_coords.size(); i++)
      this->vertex_coords.push_back(vertex_coords[i]);

    this->limit_cell.reserve(limit_cell.size());
    for (int i = 0; i < limit_cell.size(); i++)
      this->limit_cell.push_back( (short) limit_cell[i]);
  }
};



} // end namespace dolfin

#endif
