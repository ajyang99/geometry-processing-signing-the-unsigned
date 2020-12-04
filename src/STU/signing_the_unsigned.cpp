#include "STU/signing_the_unsigned.h"
#include "STU/coarse_mesh.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream>

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
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Get a coarse tet mesh, along with the rough unsigned dist estimates
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd D;
  Eigen::VectorXi I;
  Eigen::MatrixXi T;
  coarse_mesh(P, x, h, g, V, I, T);
  std::cout << "number of tets: " << T.rows() << std::endl;
  int num_tets = T.rows();
  F.resize(num_tets * 4, 3);
  for (unsigned i = 0; i < num_tets; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      for (unsigned k = 0; k < 3; ++k) {
        F(i*4+j, k) = T(i, (j+k+1)%4);
        // if ((T(i, (j+k+1)%4) < 0) || (T(i, (j+k+1)%4)>= V.rows())) {
        //   std::cout << "WARNING: T " << i << " " << (j+k+1)%4 << " " << T(i, (j+k+1)%4) << std::endl;
        // }
        // if ((F(i*4+j,k) < 0) || F(i*4+j,k) >= V.rows()) {
        //   std::cout << "WARNING: F " << i*4+j << " " << k << " " << F(i*4+j,k) << std::endl;
        // }
      }
    }
  }

  // std::cout << "AFTER" << std::endl;
  // for (unsigned i = 0; i < num_tets * 4; ++i) {
  //   for (unsigned j = 0; j < 3; ++j) {
  //     if ((F(i,j) < 0) || F(i,j) >= V.rows()) {
  //       std::cout << "WARNING: F " << i << " " << j << " " << F(i,j) << std::endl;
  //     }
  //   }
  // }
  // for (unsigned i = 0; i < num_tets; ++i) {
  //   for (unsigned j = 0; j < 4; ++j) {
  //     if ((T(i,j) < 0) || T(i,j) >= V.rows()) {
  //       std::cout << "WARNING: T " << i << " " << j << " " << T(i,j) << std::endl;
  //     }
  //   }
  // }

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  // igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}