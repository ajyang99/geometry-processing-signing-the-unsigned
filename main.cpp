#include "PSR/poisson_surface_reconstruction.h"
#include "STU/signing_the_unsigned.h"
#include "STU/graph_representation.h"
#include "STU/shoot_ray.h"
#include "STU/eps_band_vis.h"
#include <igl/list_to_matrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <random>

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
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);

  // Reconstruct mesh with PSR
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd poisx;
  Eigen::VectorXd poisg;
  poisson_surface_reconstruction(P,N,V,F,poisx, poisg);

  // Reconstruct mesh with signing the unsigned
  Eigen::MatrixXd V_stu, DG, finalV;
  Eigen::MatrixXi F_stu_coarse, F_stu_eps, T, finalF;
  Eigen::VectorXd dist, signconf, signdist;
  Eigen::VectorXi sign;
  double eps, h;
  signing_the_unsigned(
    P, V_stu, F_stu_coarse, F_stu_eps, T, eps, dist, DG,
    sign, signconf, signdist, finalV, finalF, h);

  /********** For visualizing ray shooting & sign estimation ******************/
  std::vector<std::vector<int>> v2v;
  std::vector<Eigen::MatrixXd> v2vec;
  graph_representation(V_stu, T, v2v, v2vec);

  // is_in_band
  Eigen::VectorXi is_in_band = Eigen::VectorXi::Zero(V_stu.rows());
  Eigen::VectorXi band_ixes(V_stu.rows());
  {int counter = 0;
  for (int i=0; i<V_stu.rows(); i++) {
    if (dist(i) < eps) {is_in_band(i) = 1;}
    else {band_ixes(counter) = i; counter++;}
  }
  band_ixes = band_ixes.topRows(counter);}

  // for plotting sign
  Eigen::MatrixXd signcolors = Eigen::MatrixXd::Zero(V_stu.rows(), 3);
  for (int i=0; i<V_stu.rows(); i++) {
    signcolors(i, sign(i) + 1) = 1.0;
  }
  // for plotting distance
  Eigen::MatrixXd distcolors = Eigen::MatrixXd::Zero(V_stu.rows(), 3);
  double maxdist = dist.maxCoeff();
  double mindist = dist.minCoeff();
  for (int i=0; i<V_stu.rows(); i++) {
    double r,g,b;
    igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_INFERNO, (dist(i) - mindist) / (maxdist - mindist), r,g,b);
    distcolors(i, 0) = r;
    distcolors(i, 1) = g;
    distcolors(i, 2) = b;
  }
  // for plotting signconf
  Eigen::MatrixXd signconfcolors = Eigen::MatrixXd::Zero(V_stu.rows(), 3);
  for (int i=0; i<V_stu.rows(); i++) {
    double r,g,b;
    igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_INFERNO, signconf(i), r,g,b);
    signconfcolors(i, 0) = r;
    signconfcolors(i, 1) = g;
    signconfcolors(i, 2) = b;
  }
  // for plotting signed distance
  Eigen::MatrixXd sdistcolors = Eigen::MatrixXd::Zero(V_stu.rows(), 3);
  maxdist = signdist.maxCoeff();
  mindist = signdist.minCoeff();
  std::cout<<"STU range: "<<maxdist<<", "<<mindist<<std::endl;
  for (int i=0; i<V_stu.rows(); i++) {
    double r,g,b;
    igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_INFERNO, (signdist(i) - mindist) / (maxdist - mindist), r,g,b);
    sdistcolors(i, 0) = r;
    sdistcolors(i, 1) = g;
    sdistcolors(i, 2) = b;
  }
  // for plottiing poisson
  Eigen::MatrixXd poiscolors = Eigen::MatrixXd::Zero(poisx.rows(), 3);
  maxdist = poisg.maxCoeff();
  mindist = poisg.minCoeff();
  std::cout<<"poisson range: "<<maxdist<<", "<<mindist<<std::endl;
  for (int i=0; i<poisx.rows(); i++) {
    double r,g,b;
    igl::colormap(igl::ColorMapType::COLOR_MAP_TYPE_INFERNO, (poisg(i) - mindist) / (maxdist - mindist), r,g,b);
    poiscolors(i, 0) = r;
    poiscolors(i, 1) = g;
    poiscolors(i, 2) = b;
  }

  /********************* For visualizing eps band ******************************/
  double eps_ratio = 0.5;
  double eps_max = dist.maxCoeff();
  Eigen::MatrixXi F_custom_eps;

  // init
  Eigen::MatrixXd zcolors = signcolors;
  double zplane = 0.0;
  Eigen::MatrixXd zverts = V_stu;

  // for (int i=0; i<V_stu.rows(); i++) {
  //   std::cout<<dist(i)<<", "<<sign(i)<<", "<<signconf(i)<<", "<<signdist(i)<<std::endl;
  // }

  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::opengl::glfw::Viewer viewer;
  std::cout<<R"(
  P,p      view point cloud
  N,n      view mesh reconstructed with Poisson Surface Reconstruction
  S,s      view mesh reconstructed with Signing the Unsigned
)";
  const auto set_points_original = [&]()
  {
    viewer.data().clear();
    viewer.data().set_points(P,Eigen::RowVector3d(1,1,1));
    viewer.data().add_edges(P,(P+0.01*N).eval(),Eigen::RowVector3d(1,0,0));
  };
  const auto set_points = [&]()
  {
    viewer.data().clear();
    viewer.data().set_points(V_stu,zcolors);
    // viewer.data().set_points(V_stu,Eigen::RowVector3d(1,1,1));
    // // viewer.data().add_edges(V_stu,(V_stu+0.01*DG).eval(),Eigen::RowVector3d(1,0,0));
    // int counter = 0;
    // for (int i=0; i<v2vec.size(); i++) {
    //   counter += v2vec[i].size();
    // }
    // Eigen::MatrixXd vpts(counter, 3);
    // Eigen::MatrixXd vnorm(counter, 3);
    // counter = 0;
    // for (int i=0; i<v2vec.size(); i++) {
    //   for (int j=0; j<v2vec[i].rows(); j++) {
    //     vpts.row(counter) = V_stu.row(i);
    //     vnorm.row(counter) = v2vec[i].row(j);
    //     counter++;
    //   }
    // }
    // viewer.data().add_edges(vpts,(vpts+0.01*vnorm).eval(),Eigen::RowVector3d(1,0,0));
  };
  const auto draw_z_plane = [&]()
  {
    std::cout<<"ZPLANE: "<<zplane<<std::endl;
    Eigen::MatrixXd colors, pts;
    colors.resizeLike(zcolors);
    pts.resizeLike(zverts);
    int counter = 0;
    for (int i=0; i<zverts.rows(); i++) {
      if (std::abs(zverts(i, 2) - zplane) < 0.501*h) {
        pts.row(counter) = zverts.row(i);
        colors.row(counter) = zcolors.row(i);
        counter++;
      }
    }
    viewer.data().clear();
    if (counter>0) {
      pts = pts.topRows(counter);
      colors = colors.topRows(counter);
      viewer.data().set_points(pts,colors);
    }
  };
  const auto ray_view = [&] ()
  {
    Eigen::VectorXd direc(3);
    for (int k=0; k<3; k++) {
        direc(k) = distribution(generator);
    }
    direc.normalize();
    std::cout<<"DIRECTION: "<<direc<<std::endl;
    viewer.data().clear();
    viewer.data().set_mesh(V_stu,F_stu_eps);

    // choose point
    int query_ix = band_ixes(std::rand() % band_ixes.rows());
    // // choose direction
    // Eigen::VectorXd direc = Eigen::VectorXd::Zero(3);
    // direc(0) = 1.0;
    std::vector<int> traj, should_count;
    int num_hits;
    shoot_ray(v2v, v2vec, query_ix, direc, is_in_band, DG, traj, should_count, num_hits);
    std::cout<<num_hits<<std::endl;
    Eigen::MatrixXd vec0, vec1;
    Eigen::MatrixXd colors = Eigen::MatrixXd::Ones(traj.size() - 1, 3);
    vec0.resize(traj.size() - 1, 3);
    vec1.resize(traj.size() - 1, 3);
    for (int i=0; i<traj.size()-1; i++) {
      vec0.row(i) = V_stu.row(traj[i]);
      vec1.row(i) = V_stu.row(traj[i+1]);
      if (should_count[i] == 0 && should_count[i+1] == 1)
        {colors(i, 1) = 0.0; colors(i, 2) = 0.0;}
      if (should_count[i] == 1 && should_count[i+1] == 0)
        {colors(i, 0) = 0.0; colors(i, 2) = 0.0;}
      if (should_count[i] == 1 && should_count[i+1] == 1)
        {colors(i, 0) = 0.0; colors(i, 1) = 0.0;}
    }
    viewer.data().add_edges(vec0,vec1,colors);
    Eigen::MatrixXd nextv(1, 3);
    Eigen::MatrixXd nextn(1, 3);
    nextv.row(0) = V_stu.row(query_ix);
    nextn.row(0) = direc;
    viewer.data().add_edges(nextv, nextv + 0.01*nextn, Eigen::RowVector3d(0,0,1));

    // add graph edges
    int counter = v2vec[query_ix].rows();
    Eigen::MatrixXd vpts(counter, 3);
    Eigen::MatrixXd vnorm(counter, 3);
    counter = 0;
    for (int j=0; j<v2vec[query_ix].rows(); j++) {
      vpts.row(counter) = V_stu.row(query_ix);
      vnorm.row(counter) = v2vec[query_ix].row(j);
      counter++;
    }
    viewer.data().add_edges(vpts,(vpts+0.01*vnorm).eval(),Eigen::RowVector3d(1,0,0));
  };
  set_points_original();
  viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key,int)
  {
    switch(key)
    {
      case 'P':
      case 'p':
        set_points_original();
        return true;
      case 'V':
      case 'v':
        set_points();
        return true;
      case 'N':
      case 'n':
        viewer.data().clear();
        viewer.data().set_mesh(V,F);
        return true;
      case 'S':
      case 's':
        viewer.data().clear();
        viewer.data().set_mesh(finalV,finalF);
        return true;
      case 'C':
      case 'c':
        viewer.data().clear();
        viewer.data().set_mesh(V_stu,F_stu_coarse);
        return true;
      case 'E':
      case 'e':
        viewer.data().clear();
        viewer.data().set_mesh(V_stu,F_stu_eps);
        return true;
      case 'K':
        std::cout<<"Increasing eps ratio from " << eps_ratio;
        eps_ratio = std::min(1.0, eps_ratio + 0.01);
        std::cout<<" to "<<eps_ratio<<std::endl;
        get_eps_band_faces(dist, T, eps_max*eps_ratio, F_custom_eps);
        viewer.data().clear();
        viewer.data().set_mesh(V_stu,F_custom_eps);
        return true;
      case 'k':
        std::cout<<"Decreasing eps ratio from " << eps_ratio;
        eps_ratio = std::max(0.0, eps_ratio - 0.01);
        std::cout<<" to "<<eps_ratio<<std::endl;
        get_eps_band_faces(dist, T, eps_max*eps_ratio, F_custom_eps);
        viewer.data().clear();
        viewer.data().set_mesh(V_stu,F_custom_eps);
        return true;
      case 'r':
        viewer.data().clear();
        ray_view();
        return true;

      case 'z':
        zplane += 0.02;
        draw_z_plane();
        return true;
      case 't':
        zplane -= 0.02;
        draw_z_plane();
        return true;
      case 'w':
        viewer.data().clear();
        viewer.data().set_points(P, Eigen::RowVector3d(0,0,0));
        return true;
      case 'y':
        zcolors = sdistcolors;
        zverts = V_stu;
        draw_z_plane();
        return true;
      case 'u':
        zcolors = signcolors;
        zverts = V_stu;
        draw_z_plane();
        return true;
      case 'i':
        zcolors = signconfcolors;
        zverts = V_stu;
        draw_z_plane();
        return true;
      case 'o':
        zcolors = distcolors;
        zverts = V_stu;
        draw_z_plane();
        return true;
      case 'm':
        zcolors = poiscolors;
        zverts = poisx;
        draw_z_plane();
        return true;
      case 'b':
        viewer.data().clear();
        viewer.data().set_points(V_stu, Eigen::RowVector3d(1,1,1));
        viewer.data().add_edges(V_stu,(V_stu+0.01*DG).eval(),Eigen::RowVector3d(1,0,0));
        return true;
    }
    return false;
  };
  viewer.data().point_size = 5;
  viewer.launch();

  return EXIT_SUCCESS;
}
