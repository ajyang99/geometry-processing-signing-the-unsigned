#include "STU/eps_band_select.h"
#include "STU/delaunay_3d.h"
#include <algorithm>
#include <set>
#include <vector>
#include <igl/parallel_for.h>
#include <igl/sort.h>
#include <igl/median.h>
#include <igl/volume.h>
#include <iostream>

// Find the connected components with union find
void connected_components(
    const Eigen::MatrixXi & T,
    const std::vector<std::vector<int>> & tet_buckets,
    std::vector<int> & c)
{
    std::set<std::set<int>*> disjoint_sets;
    c.clear();
    int num_buckets = tet_buckets.size();
    c.reserve(num_buckets);
    for (const std::vector<int> & I : tet_buckets) {
        for (int i = 0; i < I.size(); ++i) {
            // We are processing T.rows(I(i))
            std::set<int>* new_set = new std::set<int>();
            std::set<int> all_vals = std::set<int>();
            for (int j = 0; j < T.cols(); ++j) {
                new_set->insert(T(I[i], j));
                all_vals.insert(T(I[i], j));
            }
            // Find all disjoint sets that intersect with the new set
            std::vector<std::set<int>*> intersected_sets;
            for (std::set<int>* set : disjoint_sets) {
                for (int val : all_vals) {
                    if (set->find(val) != set->end()) {
                        intersected_sets.emplace_back(set);
                        // it is impossible another disjoint set contains val
                        all_vals.erase(val);
                        break;
                    }
                }
                if (all_vals.size() == 0) {
                    break;
                }
            }
            // Merge all intersected sets into the new set
            for (std::set<int>* set : intersected_sets) {
                new_set->insert(set->begin(), set->end());
                disjoint_sets.erase(set);
                free(set);
            }
            disjoint_sets.insert(new_set);
        }
        c.push_back(disjoint_sets.size()); // # of connected components
    }
    for (std::set<int>* set : disjoint_sets) {
        free(set);
    }
}

// Bucket tets into buckets with unsigned dist threshold defined by eps_list
// i.e. if reverise is false, tet_buckets[:i] contains all tets whose unsigned distance of all vertices
// are <= eps_list[i];
// otherwise, tet_buckets[i] contains all tets whose unsigned distance is > eps_list[buckets-i-1] (complement)
void bucket_tets(
    const Eigen::MatrixXi & T,
    const Eigen::VectorXd & D,
    const std::vector<double> & eps_list,
    const bool reverse,
    std::vector<std::vector<int>> & tet_buckets)
{
    int num_buckets = eps_list.size();
    Eigen::VectorXi bucket_inds;
    bucket_inds.resize(T.rows());
    double sign = reverse ? -1.0 : 1.0;
    igl::parallel_for(T.rows(),[&](size_t i)
    {
        Eigen::VectorXd dist;
        dist.resize(4);
        dist << D(T(i,0)), D(T(i,1)), D(T(i,2)), D(T(i,3));
        dist *= sign;
        double tet_dist = dist.maxCoeff();
        for (int j = 0; j < num_buckets; ++j) {
            // note that eps_list is in ascending order, so if we bucket in the reversee order
            // then we have to traverse eps_list revesely
            double threshold = reverse ? eps_list[num_buckets-1-j] : eps_list[j];
            threshold *= sign;
            if ((reverse && tet_dist < threshold) || (!reverse && tet_dist <= threshold)) {
                bucket_inds(i) = j;
                break;
            }
            bucket_inds(i) = -1;
        }
    },1000);
    tet_buckets = std::vector<std::vector<int>>(num_buckets, std::vector<int>());
    for (int i = 0; i < T.rows(); ++i) {
        int j = bucket_inds(i);
        if (j != -1) {
            tet_buckets[j].emplace_back(i);
        }
    }
}

// count[:i] will store the number of elements in D_P that are <= eps[i]
void input_point_count(
    const Eigen::VectorXd & D_P,
    const std::vector<double> & eps_list,
    std::vector<int> & count)
{
    int num_buckets = eps_list.size();
    count = std::vector<int>(num_buckets, 0);
    for (int i = 0; i < D_P.size(); ++i) {
        for (int j = 0; j < num_buckets; ++j) {
            if (D_P(i) <= eps_list[j]) {
                count[j]++;
                break;
            }
        }
    }
}

void eps_band_select(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & T,
    const Eigen::VectorXd & D,
    const Eigen::VectorXd & D_P,
    double & eps)
{
    // Sort the vetices by the distance value in ascending order
    int n = D.size();
    Eigen::MatrixXd D_mat, D_sorted;
    Eigen::MatrixXi I_sorted;
    D_mat.resize(n, 1);
    D_mat.col(0) = D;  // we have to turn D into a matrix to use igl::sort
    igl::sort(D, 1, true, D_sorted, I_sorted);
    Eigen::VectorXi I_sorted_vec = I_sorted.col(0);

    // Divide the vertices into 200 equal-sized chunks based on sorted distances
    int num_buckets = 200;
    std::vector<Eigen::VectorXi> v_buckets; // v_buckets[i] has the corresponding vertex inds
    v_buckets.reserve(num_buckets);
    std::vector<double> eps_list;
    eps_list.reserve(num_buckets);
    int start_idx = 0;
    int end_idx = 0;
    for (int i = 0; i < num_buckets; ++i) {
        start_idx = end_idx;
        end_idx = int(n * 1.0 * (i+1) / num_buckets);
        Eigen::VectorXi I_block = I_sorted_vec.segment(start_idx, end_idx-start_idx);
        v_buckets.push_back(I_block);
        eps_list.push_back(D_sorted(end_idx-1, 0));
    }

    // Group the tets into the buckets
    std::vector<std::vector<int>> tet_buckets;
    bucket_tets(T, D, eps_list, false, tet_buckets);
    std::vector<std::vector<int>> tet_comp_buckets;
    bucket_tets(T, D, eps_list, true, tet_comp_buckets);
    // for (int i = 0; i < num_buckets; ++i) {
    //     std::cout << "bucket_tets " << i << " " << tet_buckets[i].size() << " " << tet_comp_buckets[i].size() << std::endl;
    // }

    // Compute M(eps)
    // std::vector<int> C_eps;  // connect components
    // connected_components(T, tet_buckets, C_eps);
    // std::vector<int> H_eps;  // cavities
    // connected_components(T, tet_comp_buckets, H_eps);
    // std::reverse(H_eps.begin(), H_eps.end()); // reverse H_eps to be consistent with C_eps

    // std::vector<double> M_eps; // (2 * (C + H) - chi) / D
    // M_eps.reserve(num_buckets);
    // V_eps.resize(0, 3); // vertices
    // T_eps.resize(0, 4); // tets
    // std::set<std::pair<int, int>> E_eps; // edges
    // std::set<std::set<int>> F_eps; // faces
    // std::vector<int> point_count;
    // input_point_count(D_P, eps_list, point_count);
    // int cum_point_count = 0;
    // for (int i = 0; i < num_buckets; ++i) {
    //     // std::cout << "Processing bucket idx " << i << std::endl;
    //     // gather all vertices, edges, faces and tets of the current eps band
    //     int verts_size_old = V_eps.rows();
    //     V_eps.conservativeResize(verts_size_old + v_buckets[i].size(), 3);
    //     for (int j = 0; j < v_buckets[i].size(); ++j) {
    //         V_eps.row(verts_size_old + j) = V.row(v_buckets[i](j));
    //     }
    //     int tets_size_old = T_eps.rows();
    //     T_eps.conservativeResize(tets_size_old + tet_buckets[i].size(), 4);
    //     for (int j = 0; j < tet_buckets[i].size(); ++j) {
    //         int t = tet_buckets[i][j];
    //         T_eps.row(tets_size_old + j) = T.row(t);
    //         for (int k = 0; k < 4; ++k) {
    //             for (int l = k + 1; l < 4; ++l) {
    //                 int start_v = T(t, k);
    //                 int end_v = T(t, l);
    //                 std::pair<int, int> new_edge = start_v < end_v ? \
    //                     std::pair<int, int>(start_v, end_v) : \
    //                     std::pair<int, int>(end_v, start_v);
    //                 E_eps.insert(new_edge);
    //             }
    //             std::set<int> new_face = std::set<int>{T(t,(k+1)%4), T(t,(k+2)%4), T(t,(k+3)%4)};
    //             F_eps.insert(new_face);
    //         }
    //     }
    //     cum_point_count += point_count[i];

    //     // std::cout << T_eps.rows() << std::endl;
    //     if (T_eps.rows() == 0) {
    //         M_eps.push_back(0);
    //         continue;
    //     }

    //     // Euler characteristic
    //     int chi = V_eps.rows() - E_eps.size() + F_eps.size() - T_eps.rows();

    //     // Density
    //     Eigen::VectorXd volumes;
    //     igl::volume(V_eps, T_eps, volumes);
    //     double density = cum_point_count * 1.0 / volumes.sum();

    //     // M
    //     double m = (2.0 * (C_eps[i] + H_eps[i]) - chi) / density;
    //     std::cout << C_eps[i] << " " << H_eps[i] << " " << chi << " " << cum_point_count << " " << volumes.sum() << std::endl;
    //     std::cout << i << m << std::endl;
    //     M_eps.push_back(m);
    // }

    // int eps_idx = num_buckets * 0.5;
    // eps = eps_list[eps_idx];
    // std::cout << eps << std::endl;
    // double eps_median;
    igl::median(D, eps);
    // std::cout << eps << std::endl;
    // V_eps.resize(n, 3);
    // I_eps.resize(n);
    // T_eps.resize(T.rows(), 4);
    // int V_eps_size = 0;
    // int T_eps_size = 0;
    // for (int i = 0; i < n; ++i) {
    //     if (D(i) <= eps) {
    //         V_eps.row(V_eps_size) = V.row(i);
    //         I_eps(V_eps_size) = i;
    //         V_eps_size++;
    //     }
    // }
    // for (int i = 0; i < T.rows(); ++i) {
    //     if ((D(T(i,0)) <= eps) && (D(T(i,1)) <= eps) && (D(T(i,2)) <= eps) && (D(T(i,3)) <= eps))  {
    //         T_eps.row(T_eps_size) = T.row(i);
    //         T_eps_size++;
    //     }
    // }
    // V_eps.conservativeResize(V_eps_size, 3);
    // I_eps.conservativeResize(V_eps_size);
    // T_eps.conservativeResize(T_eps_size, 4);
}