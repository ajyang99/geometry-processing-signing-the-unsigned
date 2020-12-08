# Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets

This repository implements
[Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets](https://hal.inria.fr/inria-00502473/document)
by Mullen et al. 2010.

**Problem:** Extract a signed distance field (SDF) given only a raw pointcloud.

**Solution:** This paper leverages the fact that, given a closed surface $S \subset \mathbb{R}^3$, a query point $\mathbf{x} \in \mathbb{R}^3$ can be classified as lying inside or outside the surface by shooting a ray in any direction from the query point and counting how many times the ray intersections the surface $S$; if the ray intersects the surface an even/odd number of times, the point is outside/inside the surface. The paper uses this fact to "sign" an unsigned distance field. The steps are:

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
* `P` / `p` - coarse mesh vertices colored according to the switches listed below
* `y` - color according to predicted signed distance
* `u` - color according to predicted sign
* `i` - color according to sign confidence
* `o` - color according to unsigned distance
* `m` - color according to poisson distance
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

## Final Refinement