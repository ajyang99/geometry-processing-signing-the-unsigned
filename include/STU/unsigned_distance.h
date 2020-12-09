#ifndef UNSIGNED_DISTANCE_H
#define UNSIGNED_DISTANCE_H
#include <Eigen/Core>

// Takes input sample points P, and sampled 3D point and gives an unsigned
// distance estimate based on k-nearest-neighbor.
//
// Inputs:
//   P      #P by 3 list of input points
//   SP     #SP by 3 list of sampled points whose unsigned distance we'd like to compute
//   k      k used in knn
// Outputs:
//   D      list of size #SP of the unsigned distance
//   DG     #SP by 3 gradient of the unsigned distance
//
void unsigned_distance(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & SP,
    size_t k,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & DG);

#endif