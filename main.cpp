#include "PSR/poisson_surface_reconstruction.h"
#include "STU/signing_the_unsigned.h"
#include <igl/list_to_matrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // Load in points + normals from .pwn file
  Eigen::MatrixXd P,N;
  {
    Eigen::MatrixXd D;
    std::vector<std::vector<double> > vD;
    std::string line;
    std::fstream in;
    in.open(argc>1?argv[1]:"../data/hand.pwn");
    while(in)
    {
      std::getline(in, line);
      std::vector<double> row;
      std::stringstream stream_line(line);
      double value;
      while(stream_line >> value) row.push_back(value);
      if(!row.empty()) vD.push_back(row);
    }
    if(!igl::list_to_matrix(vD,D)) return EXIT_FAILURE;
    assert(D.cols() == 6 && "pwn file should have 6 columns");
    P = D.leftCols(3);
    N = D.rightCols(3);
  }

  // Reconstruct mesh with PSR
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  poisson_surface_reconstruction(P,N,V,F);

  // Reconstruct mesh with signing the unsigned
  Eigen::MatrixXd V_stu;
  Eigen::MatrixXi F_stu;
  signing_the_unsigned(P,V_stu,F_stu);

  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
  P,p      view point cloud
  M,m      view mesh reconstructed with Poisson Surface Reconstruction
  S,s      view mesh reconstructed with Signing the Unsigned
)";
  const auto set_points = [&]()
  {
    viewer.data().clear();
    viewer.data().set_points(P,Eigen::RowVector3d(1,1,1));
    viewer.data().add_edges(P,(P+0.01*N).eval(),Eigen::RowVector3d(1,0,0));
  };
  set_points();
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key,int)
  {
    switch(key)
    {
      case 'P':
      case 'p':
        set_points();
        return true;
      case 'M':
      case 'm':
        viewer.data().clear();
        viewer.data().set_mesh(V,F);
        return true;
      case 'S':
      case 's':
        viewer.data().clear();
        viewer.data().set_mesh(V_stu,F_stu);
        return true;
    }
    return false;
  };
  viewer.data().point_size = 2;
  viewer.launch();

  return EXIT_SUCCESS;
}

