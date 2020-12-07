#ifndef GRAPH_REPRESENTATION_H
#define GRAPH_REPRESENTATION_H
#include <Eigen/Core>

// Takes input sample points P and gives a watertight mesh.
//
// Inputs:
//   P  #P by 3 list of input points
// Outputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh triangle indinces into V
//
void graph_representation(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & T,
    std::vector<std::vector<int>> & v2v,
    std::vector<Eigen::MatrixXd> & v2vec);

#endif