#ifndef COARSE_MESH_H
#define COARSE_MESH_H
#include <Eigen/Core>

// Takes input sample points P, and sampled 3D point and gives an unsigned
// distance estimate based on k-nearest-neighbor.
//
// Inputs:
//   P      #P by 3 list of input points
//   SP     #SP by 3 list of sampled candidate points, which is a voxel grid
//   h      width of each cubic voxel grid cell
// Outputs:
//   D      list of size #SP of the unsigned distance
//   V      #V by 3 vertices of the coarse mesh (#V <= #SP)
//   I      list of size #V of the indices in SP
//   T      #T by 4 faces of the coarse tet mesh
//
void coarse_mesh(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & SP,
    const double h,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & V,
    Eigen::VectorXi & I,
    Eigen::MatrixXi & T);

#endif