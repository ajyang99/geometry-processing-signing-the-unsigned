#ifndef EPS_BAND_VIS_H
#define EPS_BAND_VIS_H
#include <Eigen/Core>

// Return all faces of a tet mesh
void get_faces(
  const Eigen::MatrixXi & T,
  Eigen::MatrixXi & F);


// Return all tets inside an eps-band
void get_eps_band(
  const Eigen::VectorXd & D,
  const Eigen::MatrixXi & T,
  const double eps,
  Eigen::MatrixXi & T_eps);


// Return faces of all tets inside an eps-band
void get_eps_band_faces(
  const Eigen::VectorXd & D,
  const Eigen::MatrixXi & T,
  const double eps,
  Eigen::MatrixXi & F_eps);

#endif