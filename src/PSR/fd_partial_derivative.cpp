#include "PSR/fd_partial_derivative.h"
#include <vector>

int flat_idx(
  const int nx, const int ny, const int nz,
  const int i, const int j, const int k)
{
  return i + j * nx + k * nx * ny;
}

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  // we assume dir \in {0, 1, 2}
  assert((dir == 0) || (dir == 1) || (dir == 2));
  
  std::vector<Eigen::Triplet<double>> triplet_list;

  if (dir == 0) {
    D.resize((nx-1)*ny*nz, nx*ny*nz);
    triplet_list.reserve((nx-1)*ny*nz*2);
    for (unsigned i = 0; i < nx - 1; ++i) {
      for (unsigned j = 0; j < ny; ++j) {
        for (unsigned k = 0; k < nz; ++k) {
          triplet_list.push_back(Eigen::Triplet<double>(
            flat_idx(nx-1, ny, nz, i, j, k), flat_idx(nx, ny, nz, i, j, k), -1
          ));
          triplet_list.push_back(Eigen::Triplet<double>(
            flat_idx(nx-1, ny, nz, i, j, k), flat_idx(nx, ny, nz, i+1, j, k), 1
          ));
        }
      }
    }
  } else if (dir == 1) {
    D.resize(nx*(ny-1)*nz, nx*ny*nz);
    triplet_list.reserve(nx*(ny-1)*nz*2);
    for (unsigned i = 0; i < nx; ++i) {
      for (unsigned j = 0; j < ny-1; ++j) {
        for (unsigned k = 0; k < nz; ++k) {
          triplet_list.push_back(Eigen::Triplet<double>(
            flat_idx(nx, ny-1, nz, i, j, k), flat_idx(nx, ny, nz, i, j, k), -1
          ));
          triplet_list.push_back(Eigen::Triplet<double>(
            flat_idx(nx, ny-1, nz, i, j, k), flat_idx(nx, ny, nz, i, j+1, k), 1
          ));
        }
      }
    }
  } else {
    D.resize(nx*ny*(nz-1), nx*ny*nz);
    triplet_list.reserve(nx*ny*(nz-1)*2);
    for (unsigned i = 0; i < nx; ++i) {
      for (unsigned j = 0; j < ny; ++j) {
        for (unsigned k = 0; k < nz-1; ++k) {
          triplet_list.push_back(Eigen::Triplet<double>(
            flat_idx(nx, ny, nz-1, i, j, k), flat_idx(nx, ny, nz, i, j, k), -1
          ));
          triplet_list.push_back(Eigen::Triplet<double>(
            flat_idx(nx, ny, nz-1, i, j, k), flat_idx(nx, ny, nz, i, j, k+1), 1
          ));
        }
      }
    }
  }

  D.setFromTriplets(triplet_list.begin(), triplet_list.end());
  D = D * (1/h);
}
