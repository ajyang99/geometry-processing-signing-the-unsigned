#ifndef EPS_BAND_SELECT_H
#define EPS_BAND_SELECT_H
#include <Eigen/Core>

// Given the coarse mesh, choose the epsilon value and return the corrresponding
// epsilon-band.
//
// Inputs:
//   V      #V by 3 vertices of the coarse tet mesh
//   T      #T by 4 facets of the coarse tet mesh
//   D      #V vector of estimated unsigned distance to the point cloud
//   D_P    #P vector of estimated unsigned distance of the sampled grid point associated to input point
// Outputs:
//   eps    chosen epsilon value
//
void eps_band_select(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & T,
    const Eigen::VectorXd & D,
    const Eigen::VectorXd & D_P,
    double & eps);

#endif