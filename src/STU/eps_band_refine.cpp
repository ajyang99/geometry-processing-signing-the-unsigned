#include "STU/eps_band_select.h"
#include <igl/knn.h>
#include <igl/octree.h>
#include <igl/parallel_for.h>
#include <igl/randperm.h>
#include <fstream>

void find_points_in_eps_band(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & D,
    const double eps,
    Eigen::MatrixXd & P_eps,
    Eigen::VectorXi & I_eps)
{
    P_eps.resizeLike(P);
    I_eps.resize(P.rows());
    int num_in_the_band = 0;
    for (int i = 0; i < P.rows(); ++i) {
        if (D(i) <= eps) {
            P_eps.row(num_in_the_band) = P.row(i);
            I_eps(num_in_the_band) = i;
            num_in_the_band++;
        }
    }
    P_eps.conservativeResize(num_in_the_band, 3);
    I_eps.conservativeResize(num_in_the_band);
}

void eps_band_refine(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & D_P,
    const Eigen::MatrixXd & V,
    const double eps,
    const size_t k,
    Eigen::VectorXd & D)
{
    // Select all input points and vertices within the epsilon band
    Eigen::MatrixXd P_eps, V_eps;
    Eigen::VectorXi tmp, I_eps;
    find_points_in_eps_band(P, D_P, eps, P_eps, tmp);
    find_points_in_eps_band(V, D, eps, V_eps, I_eps);

    // Perform KNN with only input points in the band
    // Convert P_eps into an octree as required by knn
    std::vector<std::vector<unsigned>> point_indices;
    Eigen::MatrixXi CH;
    Eigen::MatrixXd CN;
    Eigen::VectorXd W;
    igl::octree(P_eps, point_indices, CH, CN, W);

    // Perform knn for V_eps
    Eigen::MatrixXi knn_indices;
    if (P_eps.rows() < k) {
        return;
    }
    igl::knn(V_eps, P_eps, k, point_indices, CH, CN, W, knn_indices);

    std::ofstream before("dist_before.txt");
    for (int i = 0; i < I_eps.size(); ++i) {
        before << I_eps(i) << " " << D(I_eps(i)) << std::endl;
    }
    before.close();

    // Refine D
    int m = 0.5 * k;
    int beta = 0.75 * k;
    igl::parallel_for(I_eps.size(),[&](size_t i)
    {
        Eigen::RowVector3d point_of_interest = V_eps.row(i);
        // We take m random subsets of size beta and fit the best plane to 
        double best_fit_residual, best_fit_dist;
        for (unsigned j = 0; j < m; ++j) {
            // Sample beta of k neighbors
            Eigen::VectorXi rand;
            igl::randperm(k, rand);
            Eigen::MatrixXd sampled_neighbors;
            sampled_neighbors.resize(beta, 3);
            for (unsigned l = 0; l < beta; ++l) {
                sampled_neighbors.row(l) = P_eps.row(knn_indices(i, rand(l)));
            }
            // Fit a plane with the neighbors
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(
                sampled_neighbors.transpose() * sampled_neighbors);
            Eigen::VectorXd u = eigensolver.eigenvectors().col(2);
            Eigen::VectorXd v = eigensolver.eigenvectors().col(1);
            Eigen::MatrixXd on_plane = (sampled_neighbors * u) * u.transpose() + (sampled_neighbors * v) * v.transpose();
            double residual = (sampled_neighbors - on_plane).rowwise().norm().mean();
            if ((j == 0) || (best_fit_residual > residual)) {
                best_fit_residual = residual;
                // update D with point to plane distance
                Eigen::RowVector3d point_on_plane = point_of_interest.dot(u) * u + point_of_interest.dot(v) * v;
                D(I_eps(i)) = (point_of_interest - point_on_plane).norm();
            }
        }
    },1000);

    std::ofstream after("dist_after.txt");
    for (int i = 0; i < I_eps.size(); ++i) {
        after << I_eps(i) << " " << D(I_eps(i)) << std::endl;
    }
    after.close();
}