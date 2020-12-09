#ifndef GRAPH_REPRESENTATION_H
#define GRAPH_REPRESENTATION_H
#include <Eigen/Core>

void graph_representation(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & T,
    std::vector<std::vector<int>> & v2v,
    std::vector<Eigen::MatrixXd> & v2vec);

#endif