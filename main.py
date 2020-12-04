import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from fire import Fire
import numpy as np
import torch
from simple_3dviz import Spherecloud, Mesh
from simple_3dviz.window import show
from sklearn.neighbors import NearestNeighbors
from pyhull.delaunay import DelaunayTri
from simple_3dviz.utils import render
from simple_3dviz.behaviours.io import SaveFrames
from simple_3dviz.behaviours.trajectory import Circle
from simple_3dviz.behaviours.movements import CameraTrajectory
from numba import jit


def read_file(fname):
    data = []
    with open(fname, 'r') as reader:
        for line in reader:
            row = list(map(float, line.split(' ')))
            data.append(row)
    return torch.Tensor(data)[:, :3]


def get_grid(data, ndisc, buffer=0.1):
    # size in meters of the loose bounding box
    diffs = data.max(0).values - data.min(0).values
    xlo = data.min(0).values - diffs * buffer
    xhi = data.max(0).values + diffs * buffer

    # discretization in meters
    grid = torch.empty(ndisc, ndisc, ndisc, 3)
    dx = []
    for i in range(3):
        dim = [1, 1, 1]
        dim[i] = ndisc
        vals = torch.linspace(xlo[i], xhi[i], ndisc)
        dx.append((vals[1] - vals[0]).item())
        grid[:, :, :, i] = vals.view(*dim).expand(ndisc, ndisc, ndisc)
    return grid.view(-1, 3), dx


def get_trivertices(tetvertices):
    return torch.cat((
        tetvertices[:, [0, 1, 2]],
        tetvertices[:, [1, 2, 3]],
        tetvertices[:, [2, 3, 0]],
        tetvertices[:, [3, 0, 1]],
    ), 0)


# def estimate_sign(ix, points, vertices, is_in_band, R):
#     for _ in range(R):
#         traj = []
#         vec = np.random.normal(size=3)
#         vec / np.linalg.norm(vec)
        
#         traj.append(is_in_band[ix])
#         while True:


# def get_v2v(vertices):
#     """v (int) -> set(v0, v1,...) (ints)
#     """
#     print('hello', len(np.unique(vertices.flatten())))
#     print(np.max(vertices) + 1)
#     # return v2v


def sign_unsigned(fname='./data/elephant.pwn', K=5, ndisc=50,
                  debug_viz=True, R=10):
    # read in data (N x 3)
    print('reading data...')
    data = read_file(fname)

    # discretize volume (XYZ x 3) and voxel size (3)
    print('voxel grid...')
    xyzpts, dx = get_grid(data, ndisc=ndisc)

    # calculate mean distance to K nearest neighbors
    print('nearest neighbors...')
    neigh = NearestNeighbors(n_neighbors=K).fit(data.numpy())
    dist, _ = neigh.kneighbors(xyzpts.numpy(), return_distance=True)
    dist = torch.Tensor(dist).mean(1)

    # filter out points that are "far" from the point cloud
    kept = dist < np.linalg.norm(dx)
    xyzpts = xyzpts[kept]
    dist = dist[kept]

    # delaunay triangulation
    print('delaunay...')
    M = DelaunayTri(xyzpts.numpy())
    points = torch.Tensor(M.points)
    tetvertices = torch.LongTensor(M.vertices)
    trivertices = get_trivertices(tetvertices)

    # choose epsilon
    eixes = dist.sort().values[torch.linspace(0, len(dist) - 1, 50).long()]
    epsilon = eixes[int(len(eixes) * 0.5)]

    # coarse sign estimate
    # print('sign estimate...')
    # v2v = get_v2v(vertices.numpy())
    # outside_ixes = torch.arange(len(dist))[dist >= epsilon]
    # is_in_band = dist < epsilon
    # for ix in outside_ixes:
    #     estimate_sign(ix.item(), points.numpy(), vertices.numpy(), is_in_band.numpy(), R)


    ### Visualization

    def debug_viz():
        print('plotting...')
        # eixes = dist.sort().values[torch.linspace(0, len(dist) - 1, 50).long()]
        # for epi, epsilon in enumerate(eixes):
        kept = (dist[trivertices[:, 0]] < epsilon) & (dist[trivertices[:, 1]] < epsilon) & (dist[trivertices[:, 2]] < epsilon)
        # if kept.sum().item() < 1:
        #     continue
        # viewdata = torch.cat((xyzpts, data), 0)
        # print(viewdata.shape, data.shape)
        # s = Spherecloud(viewdata.numpy(), sizes=[0.005 for _ in range(viewdata.shape[0])])

        m = Mesh.from_faces(points.numpy(), trivertices.numpy()[kept], colors=np.ones((len(points), 3))*[1.0, 0.0, 0.0])
        c = Circle(center=(0, 0, 0), point=(5, 0, 0.), normal=(0, 0, 1))
        ctrj = CameraTrajectory(c, speed=0.005)
        # render(m, behaviours=[ctrj, SaveFrames(f"./frame_{epi:03}_{{:03d}}.png", every_n=5)], n_frames=100)
        show(m)

    if debug_viz:
        debug_viz()


if __name__ == '__main__':
    Fire({
        'sign_unsigned': sign_unsigned,
    })
