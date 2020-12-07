#ifndef EPS_BAND_REFINE_H
#define EPS_BAND_REFINE_H
#include <Eigen/Core>

// Refine the unsigned distance inside the epsilon band.
//
// Inputs:
//   P      #P by 3 list of the input point cloud
//   D_P    #P vector of unsigned distance of the input points
//   V      #V by 3 vertices of the coarse mesh
//   eps    epsilon value for the eps-band
//   k      k used in KNN for refining estimated unsigned distance
// Outputs:
//   D      #V vector of refined unsigned distance
//
void eps_band_refine(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & D_P,
    const Eigen::MatrixXd & V,
    const double eps,
    const size_t k,
    Eigen::VectorXd & D);

#endif