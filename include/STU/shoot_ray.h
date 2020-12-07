#ifndef SHOOT_RAY_H
#define SHOOT_RAY_H
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
//
void shoot_ray(
    const std::vector<std::vector<int>> & v2v,
    const std::vector<Eigen::MatrixXd> & v2vec,
    const int query_ix,
    const Eigen::VectorXd & direc,
    const Eigen::VectorXi & is_in_band,
    const Eigen::MatrixXd & DG,
    std::vector<int> & traj,
    std::vector<int> & should_count,
    int & num_hits);

#endif