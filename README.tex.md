# Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets

**Authors:** Joyce Yang and Jonah Philion

**Project:** This repository implements
[Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets](https://hal.inria.fr/inria-00502473/document)
by Mullen et al. 2010.

**Challenge:** Extract a signed distance field (SDF) given only a raw pointcloud.

**Proposed Solution:** This paper leverages the fact that, given a closed surface $S \subset \mathbb{R}^3$, a query point $\mathbf{x} \in \mathbb{R}^3$ can be classified as lying inside or outside the surface by shooting a ray in any direction from the query point and counting how many times the ray intersections the surface $S$; if the ray intersects the surface an even/odd number of times, the point is outside/inside the surface. The paper uses this fact to "sign" an unsigned distance field. The steps are:

	1. calculate an unsigned distance field (UDF)
	2. recover a coarse mesh (M) from the UDF
	3. shoot rays to sign the UDF relative to M
	4. smooth the signed distance

**Takeaways:** Our main critique of this paper is that it has a chicken-and-the-egg problem; ray-shooting only returns reasonable estimates of the sign if the coarse mesh M is good, but if the coarse mesh M is good then the performance benefit of estimating the sign will be small. We find the algorithm is very sensitive to its many hyper-parameters, a problem that follow-up work to this paper attempts to address ([Noise-Adaptive Shape Reconstruction from Raw Point Sets](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.679.2055&rep=rep1&type=pdf)). 

> **To get started:** Clone this repository recursively
> 
>     git clone --recursive http://github.com/ajyang99/geometry-processing-signing-the-unsigned.git
>

## Installation and Compilation

We use `CGAL` for 3D Delaunay triangulation. `CGAL` requires `gmp`, `mpfr` and `boost`. If you are a Mac user,
install these dependencies with `brew install`. If you are a Windows user,
you will only need to install `boost` as the makefiles download and install the other two
for you automatically.

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
* `V` / `v` - intermediate coarse mesh output by "Signing the Unsigned"
* `N` / `n` - final mesh output by Poisson Surface Reconstruction
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
oriented normals as input.

At a high level, for a noisy input point set $\mathbf{p} \in \mathbf{P}$ the algorithm first discretizes the bounding box... TODO

## Unsigned Distance Estimation

## Coarse Mesh Construction

## Epsilon-Band Selection

## Epsilon-Band Refinement

## Sign Estimates
We first preprocess a "graph-like" representation of the coarse mesh such that for any vertex $v$, we have the edges connected to $v$, a unit vector representing the direction of each edge, and the direction of the gradient of the unsigned distance at that vertex

    $$ \nabla d_U(x) \propto \frac{1}{K} \sum_{i=1}^K (x - x_i) $$

To "shoot rays", we sample a random unit vector `direc`, then choose the edge that most closely aligns with `direc` until we reach a vertex that was already visited. We shoot $R=15$ rays from each vertex.

We count how many times each trajectory passes through the epsilon band by counting how many times the ray transitions from a vertex outside the band to a vertex inside the band and later departs from the band. As suggested in the paper, to somewhat filter out cases where a ray "grazes" the epsilon band and comes out the same side that it entered, we only count cases where the dot product of the gradient of the unsigned distance at the entrance point has a negative dot product with the gradient of the unsigned distance at the exit point.
* A ray shot from inside the elephant's belly intersects the epsilon band 3 times and is therefore correctly marked as an "interior" point. Entrances are colored green and exits are colored red.
![](images/rayshoot.png)
* Gradient of the unsigned distance. Since it roughly aligns with surface normals, the unsigned gradient is use to filter band-ray intersections in which the ray exits the same side of the band that it entered.
![](images/udistgrad.png)
* An example where a ray "grazes" the epsilon band (black edges) and our algorithm (correctly) does not count it because the gradient of the unsigned distance at the entrance aligns too closely with the gradient at the exit. 
![](images/shallowhit.png)
* Final predicted sign. Blue corresponds to "exterior" points, green corresponds to "band" points, and red corresponds to "interior" points.
![](images/finalsign.png)
* Final "confidence" in the predicted sign. Dark corresponds to low confidence and bright corresponds to high confidence. The ears are regions of low confidence, which is intuitive because they are thin.
![](images/finalconf.png)

## Final Refinement
We now have unsigned distances for all vertices in the coarse mesh and we have an estimate of the sign of vertices outside the epsilon band (along with a measure of uncertainty in our estimate).
Before applying marching tetrahedra, the paper suggests that we first propagate sign estimates to the "band" vertices.

To so, we sort the band vertices by their unsigned distance. Starting with the vertex with greatest distance, we check all vertices that are connected to the current vertex that have been assigned a sign. If all neighbors have the same sign, we set the the sign of that vertex to the sign of its neighbors and set the confidence to the maximum of the confidence of its neighbors.

* Before and after sign propagation.
![](images/beforeprop.png)
![](images/afterprop.png)

Before feeding our SDF to marching tets, we smooth the signed distances by solving the sparse linear system
$$ (L + \alpha W) F = \alpha W \bar{\Lambda} \bar{F} $$
where $L$ is the laplacian of the coarse mesh, $\alpha \in \mathbb{R}$ is a scalar that controls how much smoothing we would like, $W$ is a diagonal matrix holding the sign confidence, $F$ is the smoothed signed distance that we solve for, $\bar{\Lambda}$ is a diagonal matrix containing the predicted sign, and $\bar{F}$ is the unsigned distance. Note that the larger alpha is, the more closely $F$ will match $\bar{\Lambda} \bar{F}$. Large $\alpha$ therefore corresponds to no smoothing. We found $\alpha=100.0$ to work best.
* Smoothed signed distance

![](images/smoothedsdist.png)

Finally, we apply marching tetrahedra to the smoothed signed distance on the coarse mesh to get a final triangle mesh.

![](images/felephant.png)
![](images/fhand.png)
![](images/fsphere.png)
