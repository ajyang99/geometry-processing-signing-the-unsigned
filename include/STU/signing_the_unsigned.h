#ifndef SIGNING_THE_UNSIGNED_H
#define SIGNING_THE_UNSIGNED_H
#include <Eigen/Core>

// Takes input sample points P and gives a watertight mesh.
//
// Inputs:
//   P  #P by 3 list of input points
// Outputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh triangle indinces into V
//
void signing_the_unsigned(
    const Eigen::MatrixXd & P,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & T,
    double & eps,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & DG,
    Eigen::VectorXi & sign,
    Eigen::VectorXd & signconf,
    Eigen::VectorXd & signdist);

#endif