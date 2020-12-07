#include "STU/eps_band_select.h"
#include "STU/delaunay_3d.h"
#include <algorithm>
#include <set>
#include <vector>
#include <igl/parallel_for.h>
#include <igl/sort.h>
#include <igl/volume.h>
#include <iostream>
#include <fstream>

// For each (cumulative) bucket of tets, find the connected components of
// the tet mesh with union find.
void connected_components(
    const Eigen::MatrixXi & T,
    const std::vector<std::vector<int>> & tet_buckets,
    std::vector<int> & c)
{
    // We will parition all referenced vertices in T into disjoint sets
    // where each set is a connected component with a unique set_id.
    // disjoint_sets[set_id] corresponds to all vertex ids in the set
    std::vector<std::vector<int>*> disjoint_sets;
    std::vector<int> set_id_per_vertex(T.maxCoeff()+1, -1);
    int num_buckets = tet_buckets.size();
    c.clear();
    c.reserve(num_buckets);
    int num_connected = 0;
    for (const std::vector<int> & I : tet_buckets) {
        for (int i = 0; i < I.size(); ++i) {
            // We are processing T.rows(I[i))
            std::set<int> intersect_set_ids;
            for (int j = 0; j < T.cols(); ++j) {
                if (set_id_per_vertex[T(I[i], j)] != -1) {
                    intersect_set_ids.insert(set_id_per_vertex[T(I[i], j)]);
                }
            }
            if (intersect_set_ids.size() == 0) {
                // The tet is not connected to any existing components,
                // so we insert a new set
                int new_set_id = disjoint_sets.size();
                std::vector<int>* new_set = new std::vector<int>(4, 0);
                for (int j = 0; j < T.cols(); ++j) {
                    (*new_set)[j] = T(I[i], j);
                    set_id_per_vertex[T(I[i], j)] = new_set_id;
                }
                disjoint_sets.emplace_back(new_set);
                num_connected++;
            } else {
                // Find the intersected set with the largest size
                int set_id = -1;
                for (const int inter_set_id : intersect_set_ids) {
                    if ((set_id == -1) || \
                        (disjoint_sets[set_id]->size() < disjoint_sets[inter_set_id]->size())) {
                        set_id = inter_set_id;
                    }
                }
                // Merge all intersected sets and the new tet
                for (const int inter_set_id : intersect_set_ids) {
                    if (set_id == inter_set_id) {
                        continue;
                    }
                    for (const int vertex : *(disjoint_sets[inter_set_id])) {
                        disjoint_sets[set_id]->emplace_back(vertex);
                        set_id_per_vertex[vertex] = set_id;
                    }
                    num_connected--;
                    free(disjoint_sets[inter_set_id]);
                    disjoint_sets[inter_set_id] = NULL;
                }
                for (int j = 0; j < T.cols(); ++j) {
                    if (set_id_per_vertex[T(I[i], j)] == set_id) {
                        continue;
                    }
                    disjoint_sets[set_id]->emplace_back(T(I[i], j));
                    set_id_per_vertex[T(I[i], j)] = set_id;
                }
            }
        }
        c.push_back(num_connected); // # of connected components
    }
    for (std::vector<int>* set : disjoint_sets) {
        if (set != NULL) {
            free(set);
        }
    }
}

// Bucket tets into buckets based on unsigned distance and threshold eps_list.
// i.e. tet_buckets[:i] contains all tets whose unsigned distance is <= eps_list[i];
// If reverse is true, tet_buckets[i] contains all tets whose unsigned distance is > eps_list[buckets-i-1].
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
    igl::parallel_for(T.rows(),[&](size_t i)
    {
        Eigen::VectorXd dist;
        dist.resize(4);
        dist << D(T(i,0)), D(T(i,1)), D(T(i,2)), D(T(i,3));
        double tet_dist = reverse? dist.minCoeff() : dist.maxCoeff();
        for (int j = 0; j < num_buckets; ++j) {
            if (reverse) {
                // note that eps_list is in ascending order, so if we bucket in the reverse order
                // then we have to traverse eps_list revesely
                double threshold = eps_list[num_buckets-1-j];
                if (tet_dist > threshold) {
                    bucket_inds(i) = j;
                    break;
                }
            } else {
                double threshold = eps_list[j];
                if (tet_dist <= threshold) {
                    bucket_inds(i) = j;
                    break;
                }
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
    std::vector<int> C_eps;  // connect components
    connected_components(T, tet_buckets, C_eps);
    std::vector<int> H_eps;  // cavities
    connected_components(T, tet_comp_buckets, H_eps);
    std::reverse(H_eps.begin(), H_eps.end()); // reverse H_eps to be consistent with C_eps

    std::vector<double> M_eps; // (2 * (C + H) - chi) / D
    M_eps.reserve(num_buckets);
    Eigen::MatrixXi T_eps;
    T_eps.resize(0, 4); // tets
    std::set<int> V_eps_in_T; // vertices referenced by T_eps
    std::set<std::pair<int, int>> E_eps; // edges
    std::set<std::set<int>> F_eps; // faces
    std::vector<int> point_count;
    input_point_count(D_P, eps_list, point_count);
    int cum_point_count = 0;
    for (int i = 0; i < num_buckets; ++i) {
        // gather all vertices, edges, faces and tets of the current eps band
        int tets_size_old = T_eps.rows();
        T_eps.conservativeResize(tets_size_old + tet_buckets[i].size(), 4);
        for (int j = 0; j < tet_buckets[i].size(); ++j) {
            int t = tet_buckets[i][j];
            T_eps.row(tets_size_old + j) = T.row(t);
            for (int k = 0; k < 4; ++k) {
                for (int l = k + 1; l < 4; ++l) {
                    int start_v = T(t, k);
                    int end_v = T(t, l);
                    V_eps_in_T.insert(start_v);
                    V_eps_in_T.insert(end_v);
                    std::pair<int, int> new_edge = start_v < end_v ? \
                        std::pair<int, int>(start_v, end_v) : \
                        std::pair<int, int>(end_v, start_v);
                    E_eps.insert(new_edge);
                }
                std::set<int> new_face = std::set<int>{T(t,(k+1)%4), T(t,(k+2)%4), T(t,(k+3)%4)};
                F_eps.insert(new_face);
            }
        }
        cum_point_count += point_count[i];

        if (T_eps.rows() == 0) {
            // edge case where the mesh is empty 
            M_eps.push_back(0);
            continue;
        }

        // Euler characteristic
        int chi = V_eps_in_T.size() - E_eps.size() + F_eps.size() - T_eps.rows();

        // Density
        Eigen::VectorXd volumes;
        igl::volume(V, T_eps, volumes);
        double density = cum_point_count * 1.0 / volumes.sum();

        // M
        double m = (2.0 * (C_eps[i] + H_eps[i]) - chi) / density;
        M_eps.push_back(m);
    }

    // Smooth M by peforming rolling average over [-smoothing_window_size_half:smoothing_window_size_half]
    int smoothing_window_size_half = 7;
    std::vector<double> M_eps_smoothed;
    for (int i = 0; i < num_buckets; ++i) {
        double s = 0;
        for (int j = i - smoothing_window_size_half; j <= i + smoothing_window_size_half; ++j) {
            s += M_eps[std::max(0, std::min(j, num_buckets-1))];
        }
        M_eps_smoothed.emplace_back(s / (2*smoothing_window_size_half+1));
    }
    // Find the local min after first local max
    int eps_idx = 0;
    while ((eps_idx < num_buckets - 1) && (M_eps_smoothed[eps_idx] <= M_eps_smoothed[eps_idx+1])) {
        eps_idx++;
    }
    while ((eps_idx < num_buckets - 1) && (M_eps_smoothed[eps_idx] >= M_eps_smoothed[eps_idx+1])) {
        eps_idx++;
    }

    eps = eps_list[eps_idx];
    std::cout << "choosing eps at " << eps_idx << "/" << num_buckets << " buckets, with value " << eps << std::endl;

    std::ofstream outfileM("M.txt");
    for (int i = 0; i < M_eps.size(); ++i) {
        outfileM << M_eps[i] << std::endl;
    }
    outfileM.close();
}