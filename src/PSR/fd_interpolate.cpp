#include "PSR/fd_interpolate.h"
#include <cmath>
#include <vector>

double interp_coeff(const unsigned offset, const double t) {
  if (offset == 0) {
    return 1 - t;
  } else {
    return t;
  }
}

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(P.rows() * 8);  // each tri-linear interpolation has 8 coefficients
  for (unsigned idx = 0; idx < P.rows(); ++idx) {
    Eigen::RowVector3i cube_corner_ind;
    unsigned x0 = int(std::floor((P(idx, 0) - corner(0)) / h));
    unsigned y0 = int(std::floor((P(idx, 1) - corner(1)) / h));
    unsigned z0 = int(std::floor((P(idx, 2) - corner(2)) / h));

    double xd = (P(idx, 0) - corner(0)) / h - x0;
    double yd = (P(idx, 1) - corner(1)) / h - y0;
    double zd = (P(idx, 2) - corner(2)) / h - z0;

    for (unsigned x_off = 0; x_off < 2; ++x_off) {
      for (unsigned y_off = 0; y_off < 2; ++y_off) {
        for (unsigned z_off = 0; z_off < 2; ++z_off) {
          triplet_list.push_back(Eigen::Triplet<double>(
            idx, (x0+x_off) + nx*((y0+y_off) + ny*(z0+z_off)),
            interp_coeff(x_off, xd) * interp_coeff(y_off, yd) * interp_coeff(z_off, zd)
          ));
        }
      }
    }
  }

  W.resize(P.rows(), nx*ny*nz);
  W.setFromTriplets(triplet_list.begin(), triplet_list.end());
}
