#include "PSR/fd_grad.h"
#include "PSR/fd_partial_derivative.h"
#include <vector>

void add_new_triplets(
  std::vector<Eigen::Triplet<double>>& triplet_list,
  const Eigen::SparseMatrix<double> & A,
  const int row_offset)
{
  for (unsigned idx = 0; idx < A.outerSize(); ++idx) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, idx); it; ++it) {
      triplet_list.push_back(Eigen::Triplet<double>(
        it.row() + row_offset, it.col(), it.value()
      ));
    }
  }
}

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  G.resize((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1), nx*ny*nz);
  
  Eigen::SparseMatrix<double> Gx;
  fd_partial_derivative(nx, ny, nz, h, 0, Gx);
  Eigen::SparseMatrix<double> Gy;
  fd_partial_derivative(nx, ny, nz, h, 1, Gy);
  Eigen::SparseMatrix<double> Gz;
  fd_partial_derivative(nx, ny, nz, h, 2, Gz);

  // concatenate respective sparse matrices from fd_partial_derivative
  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(Gx.nonZeros() + Gy.nonZeros() + Gz.nonZeros());
  add_new_triplets(triplet_list, Gx, 0);
  add_new_triplets(triplet_list, Gy, Gx.rows());
  add_new_triplets(triplet_list, Gz, Gx.rows() + Gy.rows());
  
  assert(G.rows() == Gx.rows() + Gy.rows() + Gz.rows());
  assert((G.cols() == Gx.cols()) && (G.cols() == Gy.cols()) && (G.cols() == Gz.cols()));
  G.setFromTriplets(triplet_list.begin(), triplet_list.end());
}
