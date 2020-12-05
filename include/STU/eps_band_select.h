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
//   V_eps  #V_eps by 3 vertices of the epsilon-band
//   I_eps  #V_eps indices of V_eps in V, i.e. V(I(i)) == V_eps(i)
//   T_eps  #T_eps by 4 tets of the epsilon-band
//
void eps_band_select(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & T,
    const Eigen::VectorXd & D,
    const Eigen::VectorXd & D_P,
    double & eps,
    Eigen::MatrixXd & V_eps,
    Eigen::VectorXi & I_eps,
    Eigen::MatrixXi & F_eps);

#endif