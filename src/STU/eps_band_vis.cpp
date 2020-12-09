#include "STU/eps_band_vis.h"

void get_faces(
  const Eigen::MatrixXi & T,
  Eigen::MatrixXi & F)
{
  int num_tets = T.rows();
  F.resize(num_tets * 4, 3);
  Eigen::MatrixXi face_order;
  face_order.resize(4, 3);
  face_order << 3, 2, 1,
                4, 3, 1,
                2, 4, 1,
                3, 4, 2;
  for (unsigned i = 0; i < num_tets; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      for (unsigned k = 0; k < 3; ++k) {
        F(i*4+j, k) = T(i, face_order(j,k)-1);
      }
    }
  }
}


void get_eps_band(
  const Eigen::VectorXd & D,
  const Eigen::MatrixXi & T,
  const double eps,
  Eigen::MatrixXi & T_eps)
{
  T_eps.resizeLike(T);
  int T_eps_size = 0;
  for (int i = 0; i < T.rows(); ++i) {
      if ((D(T(i,0)) <= eps) && (D(T(i,1)) <= eps) && (D(T(i,2)) <= eps) && (D(T(i,3)) <= eps))  {
          T_eps.row(T_eps_size) = T.row(i);
          T_eps_size++;
      }
  }
  T_eps.conservativeResize(T_eps_size, 4);
}

void get_eps_band_faces(
  const Eigen::VectorXd & D,
  const Eigen::MatrixXi & T,
  const double eps,
  Eigen::MatrixXi & F_eps)
{
    Eigen::MatrixXi T_eps;
    get_eps_band(D, T, eps, T_eps);
    get_faces(T_eps, F_eps);
}
