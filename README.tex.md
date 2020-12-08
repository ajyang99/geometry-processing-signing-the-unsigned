# Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets

This is an implementation of the
[Signing the Unsigned: Robust Surfact Reconstruction from Raw Pointsets](https://hal.inria.fr/inria-00502473/document)
paper by Mullen et al. 2010.

> **To get started:** Clone this repository then issue
> 
>     git clone --recursive http://github.com/ajyang99/geometry-processing-signing-the-unsigned.git
>

## Installation and Compilation

To build the project, run

    mkdir build
    cd build
    cmake -DLIBIGL_WITH_CGAL=ON -DCMAKE_BUILD_TYPE=Release ../
    make

Note that we will use `CGAL` for 3D Delaunay triangulation. Make sure you installed
the dependences (`gmp`, `mpfr`, `boost`) on your machines. If you are a Mac user,
install the dependencies with `brew install`. If you are a Windows user,
you will only need to install `boost` as the makefiles download and install the other two
for you automatically.

## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./mesh-reconstruction [path to point cloud]

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