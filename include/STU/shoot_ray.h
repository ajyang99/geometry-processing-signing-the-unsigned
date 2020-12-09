#ifndef SHOOT_RAY_H
#define SHOOT_RAY_H
#include <Eigen/Core>

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