#include "STU/unsigned_distance.h"
#include <cmath>
#include <vector>
#include <igl/knn.h>
#include <igl/parallel_for.h>
#include <igl/octree.h>
// #include <igl/slice.h>
#include <iostream>

void unsigned_distance(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & SP,
    size_t k,
    Eigen::VectorXd & D)
{
    // Convert P into an octree as required by knn
    std::vector<std::vector<unsigned>> point_indices;
    Eigen::MatrixXi CH;
    Eigen::MatrixXd CN;
    Eigen::VectorXd W;
    igl::octree(P, point_indices, CH, CN, W);

    // Perform knn for SP
    Eigen::MatrixXi knn_indices;
    igl::knn(SP, P, k, point_indices, CH, CN, W, knn_indices);

    // Compute unsigned distance based on the neighbors
    int n = SP.rows();
    D.resize(n);
    igl::parallel_for(n,[&](size_t i)
    {
        Eigen::MatrixXd neighbors;
        neighbors.resize(k, P.cols());
        for (unsigned j = 0; j < k; ++j) {
            neighbors.row(j) = P.row(knn_indices(i, j));
        }
        // igl::slice(P, knn_indices.row(i), 1, neighbors); // there seems to be buggy behavior with slice
        double squared_norm_mean = (neighbors.rowwise() - SP.row(i)).rowwise().squaredNorm().mean();
        D(i) = std::sqrt(squared_norm_mean);
    },1000);
}