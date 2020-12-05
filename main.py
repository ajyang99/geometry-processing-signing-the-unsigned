import matplotlib as mpl
mpl.use('Agg')
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
from tqdm import tqdm
import json


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


@jit(nopython=True)
def calculate_sign_changes(signs):
    if len(signs) < 2:
        return 0
    counter = 0
    for signi in range(len(signs) - 1):
        if signs[signi] == 0 and signs[signi + 1] == 1:
            counter += 1
    return counter


def estimate_sign(ix, points, v2v, is_in_band, R):
    guesses = []
    for r in range(R):
        traj = [ix]
        vec = np.random.normal(size=3)
        vec = vec / np.linalg.norm(vec)

        while True:
            pos = v2v[traj[-1]]
            pts = points[pos] - points[traj[-1]]
            pts = pts / np.linalg.norm(pts, axis=1, keepdims=True)
            dist = np.matmul(pts, vec)
            best_ix = np.argmax(dist)
            if (len(traj) > 1 and pos[best_ix] == traj[-2]) or (len(traj) > 2 and pos[best_ix] == traj[-3]):
                break
            traj.append(pos[best_ix])
        signs = is_in_band[traj]
        guesses.append(calculate_sign_changes(signs))
        # guesses.append(traj)
    return guesses


def get_v2v(tetvertices):
    """v (int) -> set(v0, v1,...) (ints)
    """
    assert(tetvertices.shape[1] == 4), tetvertices.shape
    V = np.max(tetvertices) + 1
    v2v = {v: [] for v in range(V)}
    for f in tetvertices:
        for fi in range(4):
            v2v[f[fi]].append(f[(fi + 1) % 4])
            v2v[f[fi]].append(f[(fi + 2) % 4])
            v2v[f[fi]].append(f[(fi + 3) % 4])
    assert(all((len(vs) > 0 for vs in v2v.values())))
    v2v = {k: list(set(v)) for k,v in v2v.items()}
    return v2v


def run_sign_estimate(chosen_ixes, points, v2v, is_in_band, R):
    guesses = {}
    for ix in tqdm(chosen_ixes):
        guesses[ix] = estimate_sign(ix, points, v2v, is_in_band, R)
    return guesses


def section1(fname='./data/elephant.pwn', K=5, ndisc=50):
    # read in data (N x 3)
    print('reading data...')
    data = read_file(fname)

    # discretize volume (XYZ x 3) and voxel size (3)
    print('voxel grid...')
    xyzpts, dx = get_grid(data, ndisc=ndisc)

    # calculate mean distance to K nearest neighbors
    print('calculate distance...')
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
    print('choosing epsilon...')
    eixes = dist.sort().values[torch.linspace(0, len(dist) - 1, 50).long()]
    epsilon = eixes[int(len(eixes) * 0.5)]
    is_in_band = dist < epsilon

    return points, tetvertices, trivertices, is_in_band, epsilon, dist


def sign_unsigned(fname='./data/elephant.pwn', K=5, ndisc=50,
                  debug_viz=True, R=15,
                  cache_sign=False, use_sign_cache=False):
    points, tetvertices, trivertices, is_in_band, epsilon, dist = section1(fname, K, ndisc)
    print('Chosen Epsilon:', epsilon.item())

    # coarse sign estimate
    print('sign estimate...')
    v2v = get_v2v(tetvertices.numpy())
    outside_ixes = torch.arange(len(dist))[dist >= epsilon]
    is_in_band = dist < epsilon

    if use_sign_cache:
        print('using sign cache...')
        with open('sign_cache.json', 'r') as reader:
            guesses = json.load(reader, parse_int=int)
            guesses = {int(k): v for k,v in guesses.items()}
    else:
        guesses = {}
        for ix in tqdm(outside_ixes.numpy()):
            guesses[ix] = estimate_sign(ix.item(), points.numpy(), v2v, is_in_band.numpy(), R)
    if cache_sign:
        print('caching sign...')
        guesses = {int(k): [int(v) for v in vec] for k,vec in guesses.items()}
        with open('sign_cache.json', 'w') as writer:
            json.dump(guesses, writer)

    print(guesses)

    # cross-sections of the sign plot
    alldata = []
    for v,changes in guesses.items():
        coords = points[v]
        vals = [val%2 for val in changes if val > 0]
        if len(vals) > 0:
            sign = 1 if np.mean(vals) > 0.5 else 0
            alldata.append([coords[0], coords[1], coords[2], np.mean(vals), sign])
    alldata = np.array(alldata)

    spec = np.sort(np.unique(alldata[:, 2]))
    for spei, spe in enumerate(spec):
        kept = alldata[:, 2] == spe

        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(1, 2)

        ax = plt.subplot(gs[0, 0])
        plt.scatter(alldata[kept, 0], alldata[kept, 1], c=alldata[kept, 3], vmin=0, vmax=1, cmap='winter')
        ax.set_aspect('equal')
        plt.xlim((np.min(alldata[:, 0]), np.max(alldata[:, 0])))
        plt.ylim((np.min(alldata[:, 1]), np.max(alldata[:, 1])))

        ax = plt.subplot(gs[0, 1])
        plt.scatter(alldata[kept, 0], alldata[kept, 1], c=alldata[kept, 4], vmin=0, vmax=1, cmap='winter')
        ax.set_aspect('equal')
        plt.xlim((np.min(alldata[:, 0]), np.max(alldata[:, 0])))
        plt.ylim((np.min(alldata[:, 1]), np.max(alldata[:, 1])))

        imname = f'out{spei:04}.jpg'
        print('saving', imname)
        plt.savefig(imname)
        plt.close(fig)

    # visualize coarse mesh
    kept = (dist[trivertices[:, 0]] < epsilon) & (dist[trivertices[:, 1]] < epsilon) & (dist[trivertices[:, 2]] < epsilon)
    m = Mesh.from_faces(points.numpy(), trivertices.numpy()[kept], colors=np.ones((len(points), 3))*[1.0, 0.0, 0.0])
    show(m)

    # # pcs = [Spherecloud(points.numpy(), sizes=0.005, colors=(0.0, 0.0, 0.0))]
    # pcs = []
    # for ix in guesses:
    #     info = []
    #     for pi, p in enumerate(guesses[ix]):
    #         color = torch.Tensor([[0.0, 0.0, 0.0] for _ in p])
    #         kept = is_in_band[p]
    #         sign_diff = calculate_sign_changes(kept.numpy())
    #         if sign_diff > 0:
    #             info.append(sign_diff % 2)
    #         # print(kept)
    #         color[kept, 0] = 1.0
    #         color[~kept, 1] = 1.0
    #         pcs.append(Spherecloud(points.numpy()[p], sizes=0.005, colors=color.numpy()))
    #     print(np.mean(info))

    # # visualization
    # print('plotting...')

    # # if kept.sum().item() < 1:
    # #     continue
    # # viewdata = torch.cat((xyzpts, data), 0)
    # # print(viewdata.shape, data.shape)

    # pts = np.concatenate((points.numpy()[inside], points.numpy()[outside]))
    # sizes = np.concatenate((np.array([0.01 for _ in inside]), np.array([0.005 for _ in outside])))
    # colors = np.concatenate((np.array([[1.0, 0, 0] for _ in inside]), np.array([[0.0, 1.0, 0] for _ in outside])))
    # s = Spherecloud(pts, sizes=sizes, colors=colors)

    # # coarse mesh
    # kept = (dist[trivertices[:, 0]] < epsilon) & (dist[trivertices[:, 1]] < epsilon) & (dist[trivertices[:, 2]] < epsilon)
    # m = Mesh.from_faces(points.numpy(), trivertices.numpy()[kept], colors=np.ones((len(points), 3))*[1.0, 0.0, 0.0])
    # pcs.append(Spherecloud(points.numpy(), sizes=0.005, colors=(0.0, 1.0, 0.0)))
    # pcs.append(m)

    # # outside band points
    # s = Spherecloud(points[inside], sizes=0.005, colors=(1.0, 0.0, 0.0))
    # render(s, behaviours=[SaveFrames(f"./frame0_{{:03d}}.png", every_n=1)], n_frames=1)

    # # inside band points
    # s = Spherecloud(points[outside], sizes=0.005, colors=(0.0, 1.0, 0.0))
    # render(s, behaviours=[SaveFrames(f"./frame1_{{:03d}}.png", every_n=1)], n_frames=1)

    # # c = Circle(center=(0, 0, 0), point=(5, 0, 0.), normal=(0, 0, 1))
    # # ctrj = CameraTrajectory(c, speed=0.005)
    
    # show(m)


if __name__ == '__main__':
    Fire({
        'sign_unsigned': sign_unsigned,
    })
