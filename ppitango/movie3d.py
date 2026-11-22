"""The split-screen 3D movie: ray-traced models with every contact and its distance, a camera that
flies in to the interface and around it, and the live 2D map alongside. Rendered with PyMOL."""
import json
import os
import shutil
import subprocess
import tempfile

import numpy as np
from pandappi.visualization.interaction_map import CHAIN_COLORS, _style

from . import movie
from .analysis import Scene
from .ensembles import _pdb_text
from .viewer import _bridges, _canonical

TAN_HALF_FOV = np.tan(np.radians(7))  # PyMOL field of view 14°


def _ease(u):
    u = np.clip(u, 0, 1)
    return u * u * (3 - 2 * u)


def _flight(s: float):
    """The camera at movie progress s (0..1): overview, fly in, rock to either side of the interface, fly out.
    Returns the turn about the upright axis, the tilt (radians) and how far it has flown in (0..1).

    The side-on view keeps the interface in sight; a full turn would park a whole chain in front of it."""
    zoom = _ease((s - 0.08) / 0.17) - _ease((s - 0.88) / 0.12)
    close = np.clip((s - 0.25) / 0.63, 0, 1)
    return np.radians(55) * np.sin(2 * np.pi * close), np.radians(15) * np.sin(np.pi * close), zoom


def _camera(s: float, d_far: float, d_near: float, reach: float) -> list:
    """PyMOL view at movie progress s (0..1)."""
    theta, tilt, zoom = _flight(s)
    d = d_far + (d_near - d_far) * zoom
    ry = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    rx = np.array([[1, 0, 0], [0, np.cos(tilt), -np.sin(tilt)], [0, np.sin(tilt), np.cos(tilt)]])
    return [*(ry @ rx).ravel(), 0, 0, -d, 0, 0, 0, max(d - reach, 1), d + reach, 0]


def _project(p: np.ndarray, view: list, width: int, height: int):
    """Pixel position of model point p under a PyMOL view (rotation stored column-major, camera on -z)."""
    pc = p @ np.array(view[:9]).reshape(3, 3) + np.array(view[9:12])
    scale = height / 2 / (-pc[2] * TAN_HALF_FOV)
    return width / 2 + pc[0] * scale, height / 2 - pc[1] * scale


def render(scene: Scene, output: str, seconds: float = 12.0, fps: int = 15, title: str = '',
           width: int = 940, height: int = 918, workers: int = 0) -> str:
    """Write the 3D movie to `output` (.mp4). Needs the `pymol` command (pymol-open-source)."""
    # the pymol command (conda, Homebrew), or the module of the pip wheel, which installs no command
    import importlib.util
    import sys
    if shutil.which('pymol'):
        pymol = [shutil.which('pymol')]
    elif importlib.util.find_spec('pymol'):
        pymol = [sys.executable, '-m', 'pymol']
    else:
        raise RuntimeError('PyMOL not found. Install it with: pip install "ppi-tango[movie3d]"  (or: conda install '
                           '-c conda-forge pymol-open-source)')
    ens = scene.ensemble
    a, b = ens.chains
    n = max(int(seconds * fps), 2)
    us = movie.positions(scene, n)
    coords = _canonical(scene)                                   # frames on the first axis
    bridges = _bridges(scene)
    at = lambda arr, u: movie.sample(np.moveaxis(arr, 0, -1), u, ens.periodic)
    shown = [k for k in range(len(scene.contacts)) if scene.occ[k].max() >= 0.5]
    reach = float(np.linalg.norm(coords[0], axis=1).max()) + 8
    face = coords[0][np.unique(np.concatenate([scene.atom_a[shown].ravel(), scene.atom_b[shown].ravel()]))]
    near = max(float(np.percentile(np.linalg.norm(face, axis=1), 90)) * 1.25, 14.0)

    with tempfile.TemporaryDirectory() as tmp:
        frames, models, labels = [], [], []
        for i, u in enumerate(us):
            xyz = at(coords, u)
            models.append(f'MODEL {i + 1}\n' + _pdb_text(ens.atoms, xyz).replace('END\n', 'ENDMDL\n'))
            f0 = int(round(u)) % ens.n_frames                    # contact atoms: those of the nearest frame
            view = _camera(i / (n - 1), reach * 0.8 / TAN_HALF_FOV, near / TAN_HALF_FOV, reach)
            contacts, texts, waters = [], [], []
            for k in shown:
                if movie.sample(scene.occ[k], u, ens.periodic) < 0.5:
                    continue
                kind = scene.contacts[k]['type']
                p, q = xyz[scene.atom_a[k, f0]], xyz[scene.atom_b[k, f0]]
                w = bridges[k][f0] if k in bridges else None
                if w is not None and not np.isnan(w[0]):      # a water-mediated contact goes through its water
                    waters.append(w.round(3).tolist())
                    contacts += [[*p.round(3).tolist(), *w.round(3).tolist(), _style(kind)[0]],
                                 [*w.round(3).tolist(), *q.round(3).tolist(), _style(kind)[0]]]
                else:
                    contacts.append([*p.round(3).tolist(), *q.round(3).tolist(), _style(kind)[0]])
                if kind in movie.DIRECTIONAL:  # distance labels are drawn on the picture afterwards, crisp
                    texts.append((*_project((p + q) / 2, view, width, height), f'{np.linalg.norm(p - q):.1f} Å',
                                  _style(kind)[0]))
            frames.append(dict(view=view, contacts=contacts, waters=waters))
            labels.append(texts)
        with open(f'{tmp}/traj.pdb', 'w') as f:
            f.write(''.join(models))
        json.dump(dict(pdb=f'{tmp}/traj.pdb', out=tmp, chains=[a, b], width=width, height=height,
                       colors={a: CHAIN_COLORS[0], b: CHAIN_COLORS[1]}, frames=frames,
                       interface={a: sorted({scene.number[(a, scene.contacts[k]['residue_a'])] for k in shown}),
                                  b: sorted({scene.number[(b, scene.contacts[k]['residue_b'])] for k in shown})}),
                  open(f'{tmp}/job.json', 'w'))
        parts = workers or 2  # PyMOL's ray tracer is itself multi-threaded
        print(f'  ray-tracing {n} frames with PyMOL ({parts} processes) ...')
        script = os.path.join(os.path.dirname(__file__), '_pymol_frames.py')
        procs = [subprocess.Popen([*pymol, '-cq', script, '--', f'{tmp}/job.json', str(p), str(parts)],
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) for p in range(parts)]
        logs = [p.communicate()[0] for p in procs]
        images = [f'{tmp}/f{i:05d}.png' for i in range(n)]
        if not all(map(os.path.exists, images)):
            raise RuntimeError('PyMOL did not render every frame:\n' + '\n'.join(logs)[-2000:])
        print('  compositing 3D movie ...')
        return movie.render(scene, output, seconds, fps, title, images=images, labels=labels)
