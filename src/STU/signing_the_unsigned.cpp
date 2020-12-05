#include "STU/signing_the_unsigned.h"
#include "STU/coarse_mesh.h"
#include "STU/eps_band_select.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <cmath>
#include <iostream>

void get_faces(
  const Eigen::MatrixXi & T,
  Eigen::MatrixXi & F)
{
  int num_tets = T.rows();
  F.resize(num_tets * 4, 3);
  Eigen::MatrixXi face_order;
  face_order.resize(4, 3);
  face_order << 3, 2, 1,
                4, 3, 1,
                2, 4, 1,
                3, 4, 2;
  for (unsigned i = 0; i < num_tets; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      for (unsigned k = 0; k < 3; ++k) {
        F(i*4+j, k) = T(i, face_order(j,k)-1);
      }
    }
  }
}

void signing_the_unsigned(
    const Eigen::MatrixXd & P,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, code from the Poisson Reconstruction Assignment
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Maximum extent (side length of bounding box) of points
  double max_extent =
    (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 30+2*pad samples
  double h  = max_extent/double(30+2*pad);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  std::cout << "nx " << nx << " ny " << ny << " nz " << nz << std::endl;
  std::cout << "corner " << corner(0) << " " << corner(1) << " " << corner(2) << std::endl;
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Get a coarse tet mesh, along with the rough unsigned dist estimates
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd D;
  Eigen::VectorXi I;
  Eigen::MatrixXi T;
  coarse_mesh(P, x, h, D, V, I, T);
  std::cout << "number of tets: " << T.rows() << std::endl;
  get_faces(T, F);

  ////////////////////////////////////////////////////////////////////////////
  // Choose an epsilon value for the epsilon band
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd D_cm; // D of the coarse mesh
  int n_cm = I.size();
  D_cm.resize(n_cm);
  for (unsigned i = 0; i < n_cm; ++i) {
    D_cm(i) = D(I(i));
  }

  Eigen::VectorXd D_P; // est unsigned dist of the sampled point closest to P
  D_P.resize(n);
  for (int i = 0; i < n; ++i) {
      Eigen::RowVector3d delta = P.row(i) - corner;
      int tx = int(std::round(delta(0)/h));
      int ty = int(std::round(delta(1)/h));
      int tz = int(std::round(delta(2)/h));
      D_P(i) = D(tx + nx*(ty + tz * ny));
  }
  Eigen::MatrixXd V_eps;
  Eigen::VectorXi I_eps;
  Eigen::MatrixXi T_eps;
  double eps = 0.0;
  eps_band_select(V, T, D_cm, D_P, eps, V_eps, I_eps, T_eps);
  V = V_eps;
  get_faces(T_eps, F);

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  // igl::copyleft::marching_cubes(D, x, nx, ny, nz, V, F);
}