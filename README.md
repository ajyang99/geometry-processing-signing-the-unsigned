# Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets

**Authors:** Anqi (Joyce) Yang and Jonah Philion

**Project:** This repository implements
[Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets](https://hal.inria.fr/inria-00502473/document)
by Mullen et al. 2010.

**Challenge:** Extract a signed distance field (SDF) given only a raw pointcloud.

**Proposed Solution:** This paper leverages the fact that, given a closed surface <img src="svgs/3f88470254efdef0c54f1426635096e6.svg?invert_in_darkmode" align=middle width=51.36978pt height=26.76201000000001pt/>, a query point <img src="svgs/d23cd14579e5d8e4f311ba70abe5da58.svg?invert_in_darkmode" align=middle width=48.493005000000004pt height=26.76201000000001pt/> can be classified as lying inside or outside the surface by shooting a ray in any direction from the query point and counting how many times the ray intersections the surface <img src="svgs/e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode" align=middle width=11.027445000000004pt height=22.46574pt/>; if the ray intersects the surface an even/odd number of times, the point is outside/inside the surface. The paper uses this fact to "sign" an unsigned distance field. The steps are:

	1. calculate an unsigned distance field (UDF)
	2. recover a coarse mesh (M) from the UDF
    3. use M to find an epsilon-band that best captures the surface boundary
	4. shoot rays through the epsilon-band to sign the UDF relative to M
	5. smooth the signed distance

**Takeaways:** Our main critique of this paper is that it has a chicken-and-the-egg problem; ray-shooting only returns reasonable estimates of the sign if the coarse mesh M is good, but if the coarse mesh M is good then the performance benefit of estimating the sign will be small. We find that the algorithm is very sensitive to its many hyper-parameters, a problem that follow-up work to this paper attempts to address ([Noise-Adaptive Shape Reconstruction from Raw Point Sets](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.679.2055&rep=rep1&type=pdf)). 

## Our Implementation

> **To get started:** Clone this repository recursively
> 
>     git clone --recursive http://github.com/ajyang99/geometry-processing-signing-the-unsigned.git
>

## Installation and Compilation

We use `CGAL` for 3D Delaunay triangulation. `CGAL` requires `gmp`, `mpfr` and `boost`.
If you are a Linux user, run

    sudo apt-get install libgmp3-dev libmpfr-dev libboost-all-dev
If you are a Mac user, run

    brew install gmp mpfr boost

If you are a Windows user,
you will only need to install `boost` as the makefiles below should
download and install the other two for you automatically.

To build the project, run

    mkdir build
    cd build
    cmake -DLIBIGL_WITH_CGAL=ON -DCMAKE_BUILD_TYPE=Release ../
    make

## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./mesh-reconstruction [path to point cloud]

## Visualization
#### Mesh Visualization
* `S` / `s` - mesh output by "Signing the Unsigned"
* `C` / `c` - intermediate coarse mesh output by "Signing the Unsigned"
* `E` / `e` - epsilon-band output by "Signing the Unsigned"
* `N` / `n` - mesh output by Poisson Surface Reconstruction
#### Ray Shooting
* `e` - sample a ray and visualize the result
#### Point Cloud Visualization
* `w` - input point cloud
* `P` / `p` - coarse mesh vertices colored based on which switch is selected below
* `y` - color according to predicted signed distance
* `u` - color according to predicted sign
* `i` - color according to sign confidence
* `o` - color according to unsigned distance
* `m` - color according to poisson distance
* `b` - show gradient of the unsigned distance (used during ray-shooting)
#### 2D Cut
* `r`- visualize a 2D cut and increment the 2D cut being visualized
* `t` - visualize a 2D cut and decrement the 2D cut being visualized


## Background
We implement the [Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets](https://hal.inria.fr/inria-00502473/document) by Mullen et al. 2010. The goal is to robustly reconstruct a closed, watertight mesh from a set of (potentially noisy) points. Unlike
Poisson Surface Reconstruction by Kazhdan et al. 2006, this algoithm does not require
oriented normals as input. At a high level, the algorithm estimates an unsigned
distance field, and signs it with ray shooting. Remember that the
[signed distance](https://en.wikipedia.org/wiki/Signed_distance_function) of
an arbitrary point <img src="svgs/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.977220000000004pt height=14.61206999999998pt/> to a solid volume <img src="svgs/9432d83304c1eb0dcb05f092d30a767f.svg?invert_in_darkmode" align=middle width=11.872245000000005pt height=22.46574pt/> with boundary <img src="svgs/762b53d40a0f9ca0e95dfee746e9c4c4.svg?invert_in_darkmode" align=middle width=21.512700000000006pt height=22.831379999999992pt/> is

<p align="center"><img src="svgs/3c7d8b3ae1b3fd6d9c7f4682d1ccae6c.svg?invert_in_darkmode" align=middle width=239.25659999999996pt height=49.31553pt/></p>

where <img src="svgs/6af96d9d7cc502df2449abbab9010da0.svg?invert_in_darkmode" align=middle width=18.744825000000002pt height=22.831379999999992pt/> is the unsigned distance, and <img src="svgs/556999321351626b005f887e937a9d1a.svg?invert_in_darkmode" align=middle width=17.746905000000005pt height=22.46574pt/> is the complement (i.e. outside the solid).

## Unsigned Distance Estimation
The paper first estimates the unsigned distance of an arbitrary point <img src="svgs/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.977220000000004pt height=14.61206999999998pt/> to the surface <img src="svgs/e257acd1ccbe7fcb654708f1a866bfe9.svg?invert_in_darkmode" align=middle width=11.027445000000004pt height=22.46574pt/> (which we wish to reconstruct) by evaluating the distance from <img src="svgs/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.977220000000004pt height=14.61206999999998pt/> to the input point set <img src="svgs/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode" align=middle width=12.922305000000005pt height=22.557149999999986pt/>.

To be robust to noises and outliers, the paper uses the measure defined in
[Geometric Inference for Measures based on Distance Functions](https://hal.inria.frinria-00383685/document) by Chazal et al. 2009, which first finds the top <img src="svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.137100000000004pt height=22.46574pt/> points in
<img src="svgs/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode" align=middle width=12.922305000000005pt height=22.557149999999986pt/> that are the closest to <img src="svgs/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.977220000000004pt height=14.61206999999998pt/>, and computes the unsigned distance
<p align="center"><img src="svgs/37e2190ed04353f5b59dacbb67cbac37.svg?invert_in_darkmode" align=middle width=226.95089999999996pt height=59.178735pt/></p>
where <img src="svgs/74b7d720587d24ac25c042bc258159d5.svg?invert_in_darkmode" align=middle width=48.64398pt height=24.65759999999998pt/> is the set of <img src="svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.137100000000004pt height=22.46574pt/>-nearest neighbors. The paper finds that choosing <img src="svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.137100000000004pt height=22.46574pt/> in the 12 to 30 range is sufficient. In our experiments with smaller and simpler point clouds, <img src="svgs/b0da3d297bbd326b47d166c182dc6290.svg?invert_in_darkmode" align=middle width=45.273855000000005pt height=22.46574pt/> gives better result.

## Coarse Mesh Construction
To discretize the space, we modify the adaptive sampling in the paper with fixed-grid
sampling, which we found sufficient in our experiments. Similar to the
[Poisson Surface Reconstruction assignment](https://github.com/alecjacobson/geometry-processing-mesh-reconstruction), we first define a regular 3D grid of voxels containing at least the bounding box of <img src="svgs/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode" align=middle width=12.922305000000005pt height=22.557149999999986pt/>. We then estimate the unsigned distance at each
sampled point in the grid, and only keep the point <img src="svgs/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.977220000000004pt height=14.61206999999998pt/> with
<p align="center"><img src="svgs/523748779f10020f720f1aa01902bded.svg?invert_in_darkmode" align=middle width=73.71803999999999pt height=16.438356pt/></p>
where <img src="svgs/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode" align=middle width=9.471165000000003pt height=22.831379999999992pt/> is the length of the voxel grid diagonal. The filtering step above aims to keep
only the set of sampled points that are close to the point cloud. With these
points as vertices, we can perform Delaunay triangulation in 3D to construct a coarse
tetrahedra mesh <img src="svgs/a27feed43ff187753bc9939ea2310530.svg?invert_in_darkmode" align=middle width=42.483045000000004pt height=24.65759999999998pt/>.

## <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-Band Selection
Since <img src="svgs/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode" align=middle width=9.471165000000003pt height=22.831379999999992pt/> is a loose threshold, the coarse mesh <img src="svgs/a27feed43ff187753bc9939ea2310530.svg?invert_in_darkmode" align=middle width=42.483045000000004pt height=24.65759999999998pt/> constructed above does not best reflect
the shape of the surface. To represent the shape better, we find an <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band
which is made up of all points, edges, faces and tetrahedra in the coarse mesh that have an unsigned distance <img src="svgs/ce9b6c510f151c6e44972cf13ed1093f.svg?invert_in_darkmode" align=middle width=24.024pt height=17.723969999999973pt/>. Our goal is thus to find the value of <img src="svgs/c88e6eec3577c600d80577992caabb22.svg?invert_in_darkmode" align=middle width=68.197965pt height=22.831379999999992pt/>
that best captures the surface boundary.

The gif below illustrates the <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band with different <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> values.

![](images/sec2_vis/eps_band_vis.gif)

To select the <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> value automatically, the paper uses a function
<p align="center"><img src="svgs/d99e19acc357bc1ad5ef863ff083de67.svg?invert_in_darkmode" align=middle width=200.49314999999999pt height=38.834894999999996pt/></p>
where <img src="svgs/9b325b9e31e85137d1de765f43c0f8bc.svg?invert_in_darkmode" align=middle width=12.924780000000005pt height=22.46574pt/>, <img src="svgs/7b9a0316a2fcd7f01cfd556eedf72e96.svg?invert_in_darkmode" align=middle width=14.999985000000004pt height=22.46574pt/>, <img src="svgs/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode" align=middle width=12.924780000000005pt height=22.46574pt/> are the number of components, cavities and tunnels in the <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band.
<img src="svgs/78ec2b7008296ce0561cf83393cb746d.svg?invert_in_darkmode" align=middle width=14.066250000000002pt height=22.46574pt/> is the density of input points in the band. The detailed steps are:

First, we sample 200 possible <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> values by sorting all vertices <img src="svgs/a9a3a4a202d80326bda413b5562d5cd1.svg?invert_in_darkmode" align=middle width=13.242075000000003pt height=22.46574pt/> based on the respective
unsigned distance, and splitting the vertices into 200 equal-sized buckets. The max unsigned
distance in each bucket will correspond to a potential <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> value.

Next, we bucket all input points <img src="svgs/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode" align=middle width=12.922305000000005pt height=22.557149999999986pt/> and tetrahedra <img src="svgs/2f118ee06d05f3c2d98361d9c30e38ce.svg?invert_in_darkmode" align=middle width=11.889405000000002pt height=22.46574pt/> in the coarse mesh into the 200
intervals based on the respective unsigned distance. We then compute the number of components of each <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band with union-find as we insert tetrahedra with increasing <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> thresholds. The number of cavities can be found similarily in the reverse order.
We then calculate
<img src="svgs/a75cc39d6a8fed15e3d956372092c18c.svg?invert_in_darkmode" align=middle width=191.06620499999997pt height=24.65759999999998pt/> where <img src="svgs/5212d968722a044959929d10373da783.svg?invert_in_darkmode" align=middle width=199.53235499999997pt height=24.65759999999998pt/> is
the Euler characteristic associated to the <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band. Finally, <img src="svgs/ea46831e409f911ae2b1de97d0d391d5.svg?invert_in_darkmode" align=middle width=33.52404pt height=24.65759999999998pt/> is
the number of input points inside the band divided by the volume of the band.

We then plot the function of <img src="svgs/577c43dbad72bced5ed631b967c67157.svg?invert_in_darkmode" align=middle width=37.1976pt height=24.65759999999998pt/>. The paper suggests smoothing the output slightly
and chooses the first local minimum after first local maximum based on empirical results.
In our experiments, we observed similar shapes of the <img src="svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.739810000000002pt height=22.46574pt/> plots as the paper, but found that
choosing <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> as the median of unsigned distance at all vertices gives the best result
with the full algorithm.

The image below shows the <img src="svgs/577c43dbad72bced5ed631b967c67157.svg?invert_in_darkmode" align=middle width=37.1976pt height=24.65759999999998pt/> function with respect to the bucket index for all three test point clouds.
Nearest-neighbor hyperparameter <img src="svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.137100000000004pt height=22.46574pt/>=5 here.
Blue denotes the raw values and orange denotes smoothed values.

![](images/sec2_vis/M_plots/all.png)

We also attach Figure 4 in the paper as reference.

![](images/sec2_vis/M_plots/paper.png)

The image below shows an visualization of the raw point set, the coarse bounding
mesh from Delaunay triangulation, and the <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band chosen based on the
auto-selected <img src="svgs/577c43dbad72bced5ed631b967c67157.svg?invert_in_darkmode" align=middle width=37.1976pt height=24.65759999999998pt/> values.

![](images/sec2_vis/meshes.png)

## <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-Band Refinement
Input points in the chosen <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band are most likely to be close to the surface and
not noisy. Therefore, we can refine the unsigned distances of vertices within the band using
only input points inside the band. Instead of directly using the <img src="svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.137100000000004pt height=22.46574pt/> nearest neighbors
inside the band, we now stochastically fit a plane with the neighbors. Specifically,
we take <img src="svgs/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.433210000000003pt height=14.155350000000013pt/> random subsets of size <img src="svgs/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode" align=middle width=10.165650000000005pt height=22.831379999999992pt/>. For each subset, we fit a plane with PCA and
computes the fitting residual as the mean point-to-plane distance in the subset.
We identify the best-fit plane as the one with the smallest residual, and 
set the refined unsigned distance as the distance between the plane and the query point
<img src="svgs/3b3d0d22282651ed3d6ad6428fcfa81d.svg?invert_in_darkmode" align=middle width=14.543430000000004pt height=24.65759999999998pt/>. The paper chooses <img src="svgs/38461bcc17683f8e0b68343101033567.svg?invert_in_darkmode" align=middle width=57.71799pt height=27.775769999999994pt/> and <img src="svgs/1a8b996b9e7c5c5407d4868e7ba4be80.svg?invert_in_darkmode" align=middle width=61.985549999999996pt height=27.775769999999994pt/>.
The figure belows shows the histogram of absolute distance between the ground-truth unsigned
distance (esitamted with the mesh from Poisson Surface Reconstruct) and
the unsigned distance estimated before and after the refinement on the elephant point cloud.
Note that after refinement we do not see an improvement. We suspect this might be because
refinement aims to address noise in the input point cloud problem while the test
elephant point cloud is very clean.

![](images/sec2_vis/refine/elephant.png)



## Sign Estimates
We first preprocess a "graph-like" representation of the coarse mesh such that for any vertex <img src="svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.557890000000002pt height=14.155350000000013pt/>, we have the edges connected to <img src="svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.557890000000002pt height=14.155350000000013pt/>, a unit vector representing the direction of each edge, and the direction of the gradient of the unsigned distance at that vertex

    $$ \nabla d_U(x) \propto \frac{1}{K} \sum_{i=1}^K (x - x_i) $$

To "shoot rays", we sample a random unit vector `direc`, then choose the edge that most closely aligns with `direc` until we reach a vertex that was already visited. We shoot <img src="svgs/517c96ccf36aee98c878c9c8a8c2c3b8.svg?invert_in_darkmode" align=middle width=50.96453999999999pt height=22.46574pt/> rays from each vertex.

We count how many times each trajectory passes through the epsilon band by counting how many times the ray transitions from a vertex outside the band to a vertex inside the band and later departs from the band. As suggested in the paper, to somewhat filter out cases where a ray "grazes" the epsilon band and comes out the same side that it entered, we only count cases where the dot product of the gradient of the unsigned distance at the entrance point has a negative dot product with the gradient of the unsigned distance at the exit point.
* A ray shot from inside the elephant's belly intersects the epsilon band 3 times and is therefore correctly marked as an "interior" point. Entrances are colored green and exits are colored red.
![](images/sec34_vis/rayshoot.png)
* Gradient of the unsigned distance. Since it roughly aligns with surface normals, the unsigned gradient is use to filter band-ray intersections in which the ray exits the same side of the band that it entered.
![](images/sec34_vis/udistgrad.png)
* An example where a ray "grazes" the epsilon band (black edges) and our algorithm (correctly) does not count it because the gradient of the unsigned distance at the entrance aligns too closely with the gradient at the exit. 
![](images/sec34_vis/shallowhit.png)
* Final predicted sign. Blue corresponds to "exterior" points, green corresponds to "band" points, and red corresponds to "interior" points.
![](images/sec34_vis/finalsign.png)
* Final "confidence" in the predicted sign. Dark corresponds to low confidence and bright corresponds to high confidence. The ears are regions of low confidence, which is intuitive because they are thin.
![](images/sec34_vis/finalconf.png)

## Final Refinement
We now have unsigned distances for all vertices in the coarse mesh and we have an estimate of the sign of vertices outside the epsilon band (along with a measure of uncertainty in our estimate).
Before applying marching tetrahedra, the paper suggests that we first propagate sign estimates to the "band" vertices.

To so, we sort the band vertices by their unsigned distance. Starting with the vertex with greatest distance, we check all vertices that are connected to the current vertex that have been assigned a sign. If all neighbors have the same sign, we set the the sign of that vertex to the sign of its neighbors and set the confidence to the maximum of the confidence of its neighbors.

* Before and after sign propagation.
![](images/sec34_vis/beforeprop.png)
![](images/sec34_vis/afterprop.png)

Before feeding our SDF to marching tets, we smooth the signed distances by solving the sparse linear system
<p align="center"><img src="svgs/8da3afabd53b75df12b03321e386af15.svg?invert_in_darkmode" align=middle width=159.87444pt height=17.598074999999998pt/></p>
where <img src="svgs/ddcb483302ed36a59286424aa5e0be17.svg?invert_in_darkmode" align=middle width=11.187330000000003pt height=22.46574pt/> is the laplacian of the coarse mesh, <img src="svgs/287c0853264d60279a73802455a2f550.svg?invert_in_darkmode" align=middle width=42.53980500000001pt height=22.64855999999997pt/> is a scalar that controls how much smoothing we would like, <img src="svgs/84c95f91a742c9ceb460a83f9b5090bf.svg?invert_in_darkmode" align=middle width=17.808285000000005pt height=22.46574pt/> is a diagonal matrix holding the sign confidence, <img src="svgs/b8bc815b5e9d5177af01fd4d3d3c2f10.svg?invert_in_darkmode" align=middle width=12.853995000000003pt height=22.46574pt/> is the smoothed signed distance that we solve for, <img src="svgs/cc368e320b51ae6247fc86fca93a7afd.svg?invert_in_darkmode" align=middle width=11.415690000000005pt height=26.97716999999998pt/> is a diagonal matrix containing the predicted sign, and <img src="svgs/e07348c1566e756fb2ca13951c3705a4.svg?invert_in_darkmode" align=middle width=12.853995000000003pt height=26.97716999999998pt/> is the unsigned distance. Note that the larger alpha is, the more closely <img src="svgs/b8bc815b5e9d5177af01fd4d3d3c2f10.svg?invert_in_darkmode" align=middle width=12.853995000000003pt height=22.46574pt/> will match <img src="svgs/13f02663c65a4fb8ec160721647c6bf9.svg?invert_in_darkmode" align=middle width=24.269520000000004pt height=26.97716999999998pt/>. Large <img src="svgs/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode" align=middle width=10.576500000000003pt height=14.155350000000013pt/> therefore corresponds to no smoothing. We found <img src="svgs/2bf52b113de97c54b8e8df7364c566cc.svg?invert_in_darkmode" align=middle width=69.93723pt height=21.18732pt/> to work best.
* Smoothed signed distance

![](images/sec34_vis/smoothedsdist.png)

Finally, we apply marching tetrahedra to the smoothed signed distance on the coarse mesh to get a final triangle mesh.

![](images/sec34_vis/felephant.png)
![](images/sec34_vis/fhand.png)
![](images/sec34_vis/fsphere.png)

## Discussion
Empirically we found that the results are sensitive to four hyperparameters, which are
the hyperparameters associated with discretization density, the nearest neighbor <img src="svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.137100000000004pt height=22.46574pt/>,
the number of rays used for sign esitmation, and the
smoothness hyperparametere <img src="svgs/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode" align=middle width=10.576500000000003pt height=14.155350000000013pt/>.
In addition, the automatic <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/>-band selection heuristics might not give overall best
result even though the band visualization with the chosen <img src="svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.672451500000003pt height=14.155350000000013pt/> looks reasonable.
A follow-up work to this paper attempts to address this problem.
([Noise-Adaptive Shape Reconstruction from Raw Point Sets](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.679.2055&rep=rep1&type=pdf)). 