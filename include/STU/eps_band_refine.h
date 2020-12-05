#ifndef EPS_BAND_REFINE_H
#define EPS_BAND_REFINE_H
#include <Eigen/Core>

// Refine the unsigned distance inside the epsilon band.
//
// Inputs:
//   P      #P by 3 list of the given point cloud
//   V_eps  #V by 3 vertices of the eps-band
// Outputs:
//   D      #V vector of refined unsigned distance
//
void eps_band_refine(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & V_eps,
    Eigen::VectorXd & D);

#endif