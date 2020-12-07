#include "STU/coarse_mesh.h"
#include "STU/delaunay_3d.h"
#include "STU/unsigned_distance.h"
#include <cmath>
#include <vector>

void coarse_mesh(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & SP,
    const double h,
    const size_t k,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & V,
    Eigen::VectorXi & I,
    Eigen::MatrixXi & T,
    Eigen::MatrixXd & DG)
{
    // Compute unsigned distance for all points in SP
    unsigned_distance(P, SP, k, D, DG);

    // Find sampled points whose unsigned distance < diagonal of voxel grid
    I.resize(SP.rows());  // reserve space generously
    int num_valid = 0;
    double threshold = std::sqrt(3 * std::pow(h, 2));
    for (unsigned i = 0; i < SP.rows(); ++i) {
        if (D(i) < threshold) {
            I(num_valid) = i;
            num_valid += 1;
        }
    }
    I.conservativeResize(num_valid);
    // igl::slice(SP, I, 1, V);
    V.resize(num_valid, SP.cols());
    for (unsigned i = 0; i < num_valid; ++i) {
        V.row(i) = SP.row(I(i));
    }

    // Reconstruct the coarse mesh with Delaunay triangulation
    delaunay_triangulation_3d(V, T);
}