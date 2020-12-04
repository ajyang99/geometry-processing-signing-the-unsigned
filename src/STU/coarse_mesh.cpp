#include "STU/coarse_mesh.h"
#include "STU/unsigned_distance.h"
#include <cmath>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <iostream>

void delaunay_triangulation_3d(
    const Eigen::MatrixXd & V,
    Eigen::MatrixXi & T)
{
   typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
   typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, Kernel> Vb;
   typedef CGAL::Triangulation_data_structure_3<Vb>                       Tds;
   typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                    Delaunay;
   typedef Kernel::Point_3                                                Point;
   std::vector< std::pair<Point,unsigned> > points(V.rows());
   for(int i = 0;i<V.rows();i++)
   {
     points[i] = std::make_pair(Point(V(i,0),V(i,1), V(i,2)),i);
   }
   Delaunay triangulation;
   triangulation.insert(points.begin(),points.end());
   T.resize(triangulation.number_of_finite_cells(),4);
   {
     int j = 0;
     for(Delaunay::Finite_cells_iterator cit = triangulation.finite_cells_begin();
         cit != triangulation.finite_cells_end(); ++cit) 
     {
       Delaunay::Cell_handle cell = cit;
       T(j,0) = cell->vertex(0)->info();
       T(j,1) = cell->vertex(1)->info();
       T(j,2) = cell->vertex(2)->info();
       T(j,3) = cell->vertex(3)->info();

       j++;
     }
   }
}


void coarse_mesh(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & SP,
    const double h,
    Eigen::VectorXd & D,
    Eigen::MatrixXd & V,
    Eigen::VectorXi & I,
    Eigen::MatrixXi & T)
{
    // Compute unsigned distance for all points in SP
    size_t k = 30;
    unsigned_distance(P, SP, k, D);
    std::cout << "min: " << D.minCoeff() << " max: " << D.maxCoeff() << std::endl;
    Eigen::VectorXd::Index index;
    D.minCoeff(&index);
    std::cout << "argmin: " << index << std::endl;
    D.maxCoeff(&index);
    std::cout << "argmax: " << index << std::endl;

    // Find sampled points if whose unsigned distance < diagonal of voxel grid
    I.resize(SP.rows());  // reserve space generously
    int num_valid = 0;
    double threshold = std::sqrt(3 * std::pow(h, 2));
    for (unsigned i = 0; i < SP.rows(); ++i) {
        if (D(i) < threshold) {
            I(num_valid) = i;
            num_valid += 1;
        }
    }

    std::cout << "num_valid " << num_valid << std::endl;
    std::cout << "num_all " << SP.rows() << std::endl;
    I.conservativeResize(num_valid);
    // igl::slice(SP, I, 1, V);
    V.resize(num_valid, SP.cols());
    for (unsigned i = 0; i < num_valid; ++i) {
        V.row(i) = SP.row(I(i));
    }

    // Reconstruct the coarse mesh with Delaunay triangulation
    delaunay_triangulation_3d(V, T);
}