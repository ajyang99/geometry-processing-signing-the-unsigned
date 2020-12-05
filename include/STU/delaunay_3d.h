#ifndef DELAUNAY_3D_H
#define DELAUNAY_3D_H
#include <Eigen/Core>

// Given vertices, construct the tet mesh with 3D Delaunay triangulation.
//
// Inputs:
//   V      #V by 3 vertices
// Outputs:
//   T      #T by 4 tets of the tet mesh
//
void delaunay_triangulation_3d(
    const Eigen::MatrixXd & V,
    Eigen::MatrixXi & T);

#endif