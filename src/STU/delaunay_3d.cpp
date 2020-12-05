#include "STU/delaunay_3d.h"

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

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