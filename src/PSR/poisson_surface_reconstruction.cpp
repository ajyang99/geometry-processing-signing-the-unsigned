#include "PSR/poisson_surface_reconstruction.h"
#include "PSR/fd_interpolate.h"
#include "PSR/fd_grad.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <igl/point_mesh_squared_distance.h>
#include <fstream>

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXd & poisx,
    Eigen::VectorXd & poisg)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
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
  double h  = max_extent/double(70+2*pad);
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
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Solve for g with the Poisson Surface Reconstruction constaints
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> Wx;
  Eigen::RowVector3d corner_x_stagger(corner(0) + h/2, corner(1), corner(2));
  fd_interpolate(nx-1, ny, nz, h, corner_x_stagger, P, Wx);
  Eigen::SparseMatrix<double> Wy;
  Eigen::RowVector3d corner_y_stagger(corner(0), corner(1)+h/2, corner(2));
  fd_interpolate(nx, ny-1, nz, h, corner_y_stagger, P, Wy);
  Eigen::SparseMatrix<double> Wz;
  Eigen::RowVector3d corner_z_stagger(corner(0), corner(1), corner(2)+h/2);
  fd_interpolate(nx, ny, nz-1, h, corner_z_stagger, P, Wz);

  Eigen::VectorXd v = Eigen::VectorXd::Zero((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1));
  v.head((nx-1)*ny*nz) = Wx.transpose() * N.col(0);
  v.segment((nx-1)*ny*nz, nx*(ny-1)*nz) = Wy.transpose() * N.col(1);
  v.tail(nx*ny*(nz-1)) = Wz.transpose() * N.col(2);

  Eigen::SparseMatrix<double> G;
  fd_grad(nx, ny, nz, h, G);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  solver.compute(G.transpose() * G);
  g = solver.solve(G.transpose() * v);

  Eigen::SparseMatrix<double> W;
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  Eigen::MatrixXd sigma = (Eigen::VectorXd::Ones(n).transpose() * W * g) / n;
  g = g - Eigen::VectorXd::Ones(nx*ny*nz) * sigma;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
  poisx = x;
  poisg = g;

  // Used for comparing STU unsigned dist estimates
  Eigen::VectorXd d;
  Eigen::VectorXi I;
  Eigen::MatrixXd C;
  igl::point_mesh_squared_distance(x, V, F, d, I, C);
  std::ofstream file("gt_dist.txt");
  for (int i = 0; i < d.size();  ++i) {
    file  <<  std::sqrt(d(i)) << std::endl;
  }
  file.close();
}
