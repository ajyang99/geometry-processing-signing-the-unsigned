#include <STU/shoot_ray.h>
#include <vector>
#include <iostream>
#include <set>

void shoot_ray(
    const std::vector<std::vector<int>> & v2v,
    const std::vector<Eigen::MatrixXd> & v2vec,
    const int query_ix,
    const Eigen::VectorXd & direc,
    const Eigen::VectorXi & is_in_band,
    const Eigen::MatrixXd & DG,
    std::vector<int> & traj,
    std::vector<int> & should_count,
    int & num_hits)
{
    std::set<int> temptraj;
    temptraj.insert(query_ix);
    traj.push_back(query_ix);

    // step through the mesh
    int counter = 0;
    while (true) {
        Eigen::MatrixXd dist = v2vec[traj[counter]] * direc;

        // find the "most aligned" vertex
        double best_dist = -100;
        int best_ix;
        for (int k=0; k<dist.rows(); k++) {
            if (dist(k) > best_dist) {
                best_ix = k;
                best_dist = dist(k);
            }
        }

        // check if we reached a boundary (e.g. we are looping)
        if (temptraj.count(v2v[traj[counter]][best_ix]))
            break;
        traj.push_back(v2v[traj[counter]][best_ix]);
        temptraj.insert(v2v[traj[counter]][best_ix]);
        // initialize all vertices as not hitting anything
        should_count.push_back(0);
        counter++;
    }

    num_hits = 0;
    if (traj.size() < 2) {
        return;
    }
    // don't feed this function vertices that are inside the band!
    if (is_in_band[traj[0]]) {
        std::cout<<"is in band!!!"<<std::endl;
        return;
    }

    int enter_ix;
    for (int i=0; i<traj.size() - 1; i++) {
        // enter band
        if (!is_in_band[traj[i]] && is_in_band[traj[i+1]]) {
            enter_ix = i+1;
        }
        // leaving band
        else if (is_in_band[traj[i]] && !is_in_band[traj[i+1]]) {
            // only count deep hits
            if (DG.row(traj[enter_ix-1]).dot(DG.row(traj[i+1])) > 0) {
                continue;
            }
            // found a region where we were inside the band
            num_hits++;
            for (int j=enter_ix; j<i+1; j++) {
                should_count[j] = 1;
            }
        }
    }

}