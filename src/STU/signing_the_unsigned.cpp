#include "STU/signing_the_unsigned.h"
#include "STU/unsigned_distance.h"
#include "STU/coarse_mesh.h"
#include "STU/eps_band_select.h"
#include "STU/shoot_ray.h"
#include "STU/graph_representation.h"
#include "STU/eps_band_refine.h"
#include <igl/doublearea.h>
#include <igl/parallel_for.h>
#include <igl/median.h>
#include <igl/cotmatrix.h>
#include <igl/marching_tets.h>
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <functional>


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
    Eigen::MatrixXi & F_coarse,
    Eigen::MatrixXi & F_eps,
    Eigen::MatrixXi & T,
    double & eps,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & DG,
    Eigen::VectorXi & sign,
    Eigen::VectorXd & signconf,
    Eigen::VectorXd & signdist,
    Eigen::MatrixXd & finalV,
    Eigen::MatrixXi & finalF,
    double & h)
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
  h = max_extent/double(40+2*pad);
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
  const size_t nearest_neighbor_k = 5;
  Eigen::MatrixXd DG_x;
  coarse_mesh(P, x, h, nearest_neighbor_k, D_x, V, I, T, DG_x);
  get_faces(T, F_coarse);

  // used for eps refinement plotting
  // std::ofstream index_file("vertex_indexing.txt");
  // for (int i = 0; i < I.size(); ++i) {
  //   index_file << I(i) << std::endl;
  // }
  // index_file.close();

  ////////////////////////////////////////////////////////////////////////////
  // Choose an epsilon value for the epsilon band, and refine unsigned dist
  // of sampled points inside the band.
  ////////////////////////////////////////////////////////////////////////////
  D.resizeLike(I);
  DG.resize(I.rows(), 3);
  for (unsigned i = 0; i < I.size(); ++i) {
    D(i) = D_x(I(i));
    DG.row(i) = DG_x.row(I(i));
  }

  Eigen::VectorXd D_P; // est unsigned dist of the input point cloud
  Eigen::MatrixXd DGnew; // new gradient of the distance (TODO use this)
  unsigned_distance(P, P, nearest_neighbor_k, D_P, DGnew);
  eps_band_select(V, T, D, D_P, eps);
  // igl::median(D, eps);

  eps_band_refine(P, D_P, V, eps, nearest_neighbor_k, D);
  Eigen::MatrixXi T_eps;  // to visulize eps band with refined distance
  get_eps_band(D, T, eps, T_eps);
  get_faces(T_eps, F_eps);

  ////////////////////////////////////////////////////////////////////////////
  // Sign the distances D of tet mesh (V, T) with eps
  ////////////////////////////////////////////////////////////////////////////
  // graph rep of the mesh
  std::vector<std::vector<int>> v2v;
  std::vector<Eigen::MatrixXd> v2vec;
  graph_representation(V, T, v2v, v2vec);
  // is_in_band
  Eigen::VectorXi is_in_band = Eigen::VectorXi::Zero(V.rows());
  for (int i=0; i<V.rows(); i++) {
    if (D(i) < eps)
      is_in_band(i) = 1;
  }

  // prepare to shoot rays
  int R = 15;  // number of rays to shoot
  sign.resize(V.rows());  // predicted sign for each vertex
  signconf.resize(V.rows());  // "confidence" in predicted sign
  // rng for choosing directions
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);

  igl::parallel_for(V.rows(),[&](size_t i) {
    // if inside the band, set sign to 0 for now
    if (is_in_band(i)) {
      sign(i) = 0;
      signconf(i) = 0.0;
    }
    else {
      std::vector<int> hit_collect;
      // INFINITE LOOP POSSIBILITY!
      while (hit_collect.size() < R) {
        // choose a random direction
        // choose direction
        Eigen::VectorXd direc(3);
        for (int k=0; k<3; k++) {
            direc(k) = distribution(generator);
        }
        direc.normalize();
        std::vector<int> traj, should_count;  // for debugging
        int num_hits;  // parity of this determines the sign
        shoot_ray(v2v, v2vec, i, direc, is_in_band, DG, traj, should_count, num_hits);
        // // only count the ray if we hit the band at least once
        // if (num_hits > 0)
        hit_collect.push_back(num_hits % 2);
      }
      // count evens and odds
      int evens = 0;
      int odds = 0;
      for (int j=0; j<hit_collect.size(); j++) {
        if (hit_collect[j] % 2 == 0) {evens++;}
        else {odds++;}
      }
      // calculate sign and confidence
      sign(i) = (evens > odds) ? 1 : -1;
      signconf(i) = 2*std::max(evens, odds) / ((float) R) - 1;
    }
  },1000);

  // solve for sign inside the band
  // - loop through vertices in the band starting with the farthest distance
  // - if neighbors that have a sign all agree about the sign, update that sign (set signconf to be max signconf of neighbors)
  {
  std::vector<std::pair<int,double>> bandvs;
  for (int i=0; i<V.rows(); i++) {
    if (is_in_band(i))
      bandvs.push_back(std::make_pair(i,D(i)));
  }
    std::sort(bandvs.begin(), bandvs.end(), [](std::pair<int,double> a, std::pair<int,double> b) {return a.second > b.second; });
    for (int i=0; i<bandvs.size(); i++) {
      int current_ix = bandvs[i].first;
      int current_sign = 0;
      double current_conf = 0.0;
      // check if all the vertices connected to this vertex with a sign agree on the sign
      for (int j=0; j<v2v[current_ix].size(); j++) {
        if (current_sign == 0) {
          current_sign = sign(v2v[current_ix][j]);}
        else if (sign(v2v[current_ix][j]) != 0 && sign(v2v[current_ix][j]) != current_sign) {
          current_sign = 0; break;}
        if (sign(v2v[current_ix][j]) != 0) {
          current_conf = std::max(current_conf, signconf(v2v[current_ix][j]));
        }
      }
      if (current_sign != 0) {
        signconf(current_ix) = current_conf;
        sign(current_ix) = current_sign;
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Solve for the final signed distance with smoothing constraints
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd tgt(V.rows());
  Eigen::SparseMatrix<double> Amat(V.rows(), V.rows());
  Eigen::SparseMatrix<double> L;
  double alpha = 100.0;
  for (int i=0; i<V.rows(); i++) {
    tgt(i) = alpha * signconf(i) * sign(i) * D(i);
  }
  // remove faces with 0 area
  Eigen::VectorXd areas;
  Eigen::MatrixXi faces;
  get_faces(T, faces);
  igl::doublearea(V, faces, areas);
  Eigen::MatrixXi correctT;
  correctT.resizeLike(faces);
  int counter = 0;
  for (int i=0; i<faces.rows(); i++) {
    if (areas(i) > 1e-8) {
      correctT.row(counter) = faces.row(i);
      counter++;
    }
  }
  correctT = correctT.topRows(counter);
  igl::cotmatrix(V, correctT, L);

  // weight matrix
  typedef Eigen::Triplet<double> Trip;
  std::vector<Trip> tripletList;
  for(int i=0; i<V.rows(); i++) {
    tripletList.push_back(Trip(i,i,alpha*signconf(i)));
  }
  Amat.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
  // std::cout<<L.rows()<<", "<<Amat.rows()<<", "<<L.cols()<<", "<<Amat.cols()<<std::endl;
  // std::cout<<Amat.coeffs().minCoeff()<<", "<<Amat.coeffs().maxCoeff()<<std::endl;
  // std::cout<<"L coeffs: "<<L.coeffs().minCoeff()<<", "<<L.coeffs().maxCoeff()<<std::endl;
  solver.compute(-L + Amat);
  signdist = solver.solve(tgt);
  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////

  // find sigma
  double sigma;
  {
  Eigen::VectorXd inbanddist(D.rows());
  int counter = 0;
  for (int i=0; i<D.rows(); i++) {
    if (is_in_band(i))
      inbanddist(counter) = signdist(i); counter++;
  }
  inbanddist = inbanddist.topRows(counter);
  igl::median(inbanddist, sigma);
  }

  std::cout << "marching cubes sigma is " << sigma << std::endl;
  Eigen::VectorXi is_in_coarse_mesh = Eigen::VectorXi::Zero(x.rows());
  for (int i = 0; i < I.size(); ++i) {
    is_in_coarse_mesh(I(i)) = 1;
  }
  Eigen::VectorXd g(x.rows());
  for (int i = 0; i < x.rows(); ++i) {
    if (is_in_coarse_mesh(i) == 0) {
      g(i) = 1;
    }
  }
  for (int i = 0; i < I.size(); ++i) {
    g(I(i)) = signdist(i) - sigma;
  }
  std::cout<<"THIS RANGE: "<<g.maxCoeff()<<", "<<g.minCoeff()<<std::endl;
  igl::marching_tets(V, T, signdist, sigma, finalV, finalF);
}