#include <iostream>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/la/IndexMap.h>
#include <dolfin/la/PETScMatrix.h>
#include <pybind11/pybind11.h>


using dolfin::IndexMap;
using po = dolfin::PETScObject;


std::shared_ptr<dolfin::PETScMatrix>
create_block_matrix(const dolfin::FunctionSpace & V,
                    std::size_t nnz)
{
    const auto & mesh = V.mesh();
    const auto & dm = V.dofmap();
    const auto & im = dm->index_map();
    const auto comm = mesh->mpi_comm();
    PetscErrorCode ierr;

    std::size_t n = im->size(IndexMap::MapSize::OWNED);
    std::size_t N = im->size(IndexMap::MapSize::GLOBAL);

    Mat A;

    ierr = MatCreate(comm, &A);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatCreate");

    // Set matrix size and number of non zeros

    ierr = MatSetSizes(A, n, n, N, N);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatSetSizes");

    ierr = MatSetType(A, MATAIJ);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatSetType");

    ierr = MatSeqAIJSetPreallocation(A, nnz, NULL);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatSeqAIJSetPreallocation");

    ierr = MatMPIAIJSetPreallocation(A, nnz, NULL, 0, NULL);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatMPIAIJSetPreallocation");

    // Set up local to global mapping

    std::size_t bs = 1;
    std::vector<std::size_t> ltg;
    dm->tabulate_local_to_global_dofs(ltg);

    ISLocalToGlobalMapping lgmap;
    ierr = ISLocalToGlobalMappingCreate(comm, bs, ltg.size(),
                                        (const int*) ltg.data(),
                                        PETSC_COPY_VALUES, &lgmap);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "ISLocalToGlobalMappingCreate");

    ierr = MatSetLocalToGlobalMapping(A, lgmap, lgmap);
    if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatSetLocalToGlobalMapping");

    return std::make_shared<dolfin::PETScMatrix>(A);
}


std::shared_ptr<dolfin::PETScMatrix>
matmul(const dolfin::PETScMatrix & A, const dolfin::PETScMatrix & B)
{
  PetscErrorCode ierr;
  PetscReal rval = 2;
  Mat Cpetsc;

  ierr = MatMatMult(A.mat(), B.mat(), MAT_INITIAL_MATRIX, rval, &Cpetsc);
  if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatMatMult");

  auto C = std::make_shared<dolfin::PETScMatrix>(Cpetsc);
  C->apply("insert");

  return C;
}


void matmul_reuse(const dolfin::PETScMatrix & A, const dolfin::PETScMatrix & B,
            dolfin::PETScMatrix & C)
{
  PetscErrorCode ierr;
  PetscReal rval = 2;
  Mat Cpetsc = C.mat();
  ierr = MatMatMult(A.mat(), B.mat(), MAT_REUSE_MATRIX, rval, &Cpetsc);
  if (ierr != 0) po::petsc_error(ierr, __FILE__, "MatMatMult");

  C.apply("insert");
}


PYBIND11_MODULE(SIGNATURE, m)
{
  m.def("create_block_matrix", &create_block_matrix);
  m.def("matmul", &matmul);
  m.def("matmul_reuse", &matmul_reuse);
}
