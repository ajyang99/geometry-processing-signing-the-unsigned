#ifndef COARSE_MESH_H
#define COARSE_MESH_H
#include <Eigen/Core>

// Given a point set P and sampled voxel grid SP with voxel cube of size h,
// construct a coarse tet mesh whose vertices are from the voxel grid with
// Delaunay triangulation.
//
// Inputs:
//   P      #P by 3 list of input points
//   SP     #SP by 3 list of sampled candidate points, which is a voxel grid
//   h      width of each cubic voxel grid cell
//   k      K used for K nearest neighbor when estimating unsigned distance
// Outputs:
//   D      list of size #SP of the unsigned distance
//   V      #V by 3 vertices of the coarse mesh (#V <= #SP)
//   I      list of size #V of the indices in SP
//   T      #T by 4 facets of the coarse tet mesh
//
void coarse_mesh(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & SP,
    const double h,
    const size_t k,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & V,
    Eigen::VectorXi & I,
    Eigen::MatrixXi & T,
    Eigen::MatrixXd & DG);

#endif