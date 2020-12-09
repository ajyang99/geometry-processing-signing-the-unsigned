#ifndef SIGNING_THE_UNSIGNED_H
#define SIGNING_THE_UNSIGNED_H
#include <Eigen/Core>

// Takes input sample points P and gives a watertight mesh.
//
// Inputs:
//   P          #P by 3 list of input points
// Outputs:
//   V          #V by 3 list of coarse mesh vertex positions
//   F_coarse   #F_coarse by 3 list of coarse mesh triangle indices into V
//   F_eps      #F_eps by 3 list of eps band mesh triangle indices into V
//   T          #T by 4 list of coarse mesh tet indices into V
//   eps        eps value chosen for eps-band
//   D          #V list of unsigned distance estimates before final refinement
//   DG         #V by 3 gradients of unsigned distance estimates before final refinement
//   sign       #V list of sign estimated from ray shooting
//   signconf   #V list of sign confidence estimated from ray shooting
//   signdist   #V list of signed distance after final smoothing
//   finalV     #finalV by 3 list of final mesh vertex positions
//   finalF     #finalF by 3 list of final mesh triangle indices into finalV
//
void signing_the_unsigned(
    const Eigen::MatrixXd & P,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F_coarse,
    Eigen::MatrixXi & F_eps,
    Eigen::MatrixXi & T,
    double & eps,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & DG,
    Eigen::VectorXi & sign,
    Eigen::VectorXd & signconf,
    Eigen::VectorXd & signdist,
    Eigen::MatrixXd & finalV,
    Eigen::MatrixXi & finalF,
    double & h);

#endif