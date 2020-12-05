#include "STU/signing_the_unsigned.h"
#include "STU/coarse_mesh.h"
#include "STU/eps_band_select.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/parallel_for.h>
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
  // face_order << 1, 2, 3,
  //               2, 3, 4,
  //               3, 4, 1,
  //               4, 1, 2;
  for (unsigned i = 0; i < num_tets; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      for (unsigned k = 0; k < 3; ++k) {
        F(i*4+j, k) = T(i, face_order(j,k)-1);
      }
    }
  }
}

void get_eps_band(
  const Eigen::VectorXd & D,
  const Eigen::MatrixXi & T,
  const double eps,
  Eigen::MatrixXi & T_eps)
{
  T_eps.resizeLike(T);
  int T_eps_size = 0;
  for (int i = 0; i < T.rows(); ++i) {
      if ((D(T(i,0)) <= eps) && (D(T(i,1)) <= eps) && (D(T(i,2)) <= eps) && (D(T(i,3)) <= eps))  {
          T_eps.row(T_eps_size) = T.row(i);
          T_eps_size++;
      }
  }
  T_eps.conservativeResize(T_eps_size, 4);
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
  Eigen::VectorXd D_x; // unsigned distance of all sampled grid points x
  Eigen::VectorXi I; // x.row(I(i)) = V.row(i)
  Eigen::MatrixXi T; // (V, T) is the coarse mesh
  coarse_mesh(P, x, h, D_x, V, I, T);
  get_faces(T, F);

  ////////////////////////////////////////////////////////////////////////////
  // Choose an epsilon value for the epsilon band
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd D; // D of the coarse mesh; #D == #V
  D.resizeLike(I);
  for (unsigned i = 0; i < I.size(); ++i) {
    D(i) = D_x(I(i));
  }

  Eigen::VectorXd D_P; // est unsigned dist of the sampled point closest to P
  D_P.resize(n);
  igl::parallel_for(n,[&](size_t i) {
      Eigen::RowVector3d delta = P.row(i) - corner;
      int tx = int(std::round(delta(0)/h));
      int ty = int(std::round(delta(1)/h));
      int tz = int(std::round(delta(2)/h));
      D_P(i) = D(tx + nx*(ty + tz * ny));
  },1000);
  double eps = 0.0;
  eps_band_select(V, T, D, D_P, eps);
  Eigen::MatrixXi T_eps;
  // to visulize eps band
  get_eps_band(D, T, eps, T_eps);
  get_faces(T_eps, F);

  ////////////////////////////////////////////////////////////////////////////
  // Sign the distances D of tet mesh (V, T) with eps
  ////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  // igl::copyleft::marching_cubes(D, x, nx, ny, nz, V, F);
}