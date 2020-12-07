#include <STU/graph_representation.h>
#include <vector>
#include <iostream>
#include <set>


void graph_representation(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & T,
    std::vector<std::vector<int>> & v2v,
    std::vector<Eigen::MatrixXd> & v2vec) {

    std::vector<std::set<int>> v2set;
    for (int i=0; i<V.rows(); i++) {
        v2set.push_back(std::set<int>());
        v2v.push_back(std::vector<int>());
    }
    for (int i=0; i<T.rows(); i++) {
        for (int j=0; j<4; j++) {
            v2set[T(i, j)].insert(T(i, (j+1)%4));
            v2set[T(i, j)].insert(T(i, (j+2)%4));
            v2set[T(i, j)].insert(T(i, (j+3)%4));
        }
    }

    for (int i=0; i<V.rows(); i++) {
        for (int ix : v2set[i]) {
            v2v[i].push_back(ix);
        }
    }

    for (int i=0; i<V.rows(); i++) {
        Eigen::MatrixXd vec(v2v[i].size(), 3);
        int counter = 0;
        for (int ix : v2v[i]) {
            vec.row(counter) = V.row(ix);
            counter++;
        }
        vec = vec.rowwise() - V.row(i);
        vec.rowwise().normalize();
        v2vec.push_back(vec);
    }

}
