"""The interactive 3D viewer: one self-contained HTML file (3Dmol.js inlined, works offline).

The real models move through the frames with every contact drawn between its two atoms and its
distance measured live, next to the same 2D map as the movie. Embeds in slides as an <iframe>.
"""
import base64
import json
from importlib import resources

import numpy as np
from pandappi.visualization.interaction_map import BURIED_MIN, CLASS_COLOR, RESIDUE_CLASS, TYPE_STYLE, _style

from .analysis import Scene
from .movie import DIRECTIONAL

DASH = {'-': '', '--': '6 4', ':': '1.5 4', '-.': '7 3 1.5 3'}


def _view_frame(scene: Scene):
    """Centre and rotation that put chain a left, chain b right, and the interface edge-on."""
    ens = scene.ensemble
    a, b = ens.chains
    xyz = scene.points[0]
    idx = np.unique(np.concatenate([scene.atom_a.ravel(), scene.atom_b.ravel()])) if len(scene.contacts) \
        else np.arange(len(xyz))
    face = xyz[idx]
    in_a = (scene.point_chain[idx] == a)
    centre = (face[in_a].mean(0) + face[~in_a].mean(0)) / 2
    ex = face[~in_a].mean(0) - face[in_a].mean(0)
    ex /= np.linalg.norm(ex)
    flat = (face - centre) - np.outer((face - centre) @ ex, ex)
    ey = np.linalg.svd(flat, full_matrices=False)[2][0]      # longest direction within the interface: vertical
    ey -= ex * (ey @ ex)
    ey /= np.linalg.norm(ey)
    return centre, np.stack([ex, ey, np.cross(ex, ey)])


def _canonical(scene: Scene) -> np.ndarray:
    """Coordinates turned so chain a sits left, chain b right, and the interface is seen edge-on: the atoms, then
    the centres of the rings and amide groups contacts end on (`Scene.points`)."""
    centre, rot = _view_frame(scene)
    return (scene.points - centre) @ rot.T


def _bridges(scene: Scene) -> dict:
    """Bridging waters of the water-mediated contacts, in the same frame of view: contact -> (frames, 3), NaN
    where the contact is not made."""
    centre, rot = _view_frame(scene)
    return {k: (pos - centre) @ rot.T for k, pos in scene.bridge.items()}


def write(scene: Scene, output: str, title: str = '') -> str:
    """Write the viewer for `scene` to `output` (.html)."""
    ens = scene.ensemble
    a, b = ens.chains
    atoms = ens.atoms
    coords = np.clip(np.round(_canonical(scene) * 100), -32768, 32767).astype('<i2')
    shown = [k for k in range(len(scene.contacts)) if scene.occ[k].max() >= 0.5]
    drawn = set(scene.drawn)
    fill = lambda res: CLASS_COLOR[RESIDUE_CLASS.get(res.split()[0], 'Other')]
    r2 = lambda v: [round(float(x), 2) for x in v]
    bridges = _bridges(scene)
    water = lambda k: [None if np.isnan(p[0]) else r2(p) for p in bridges[k]] if k in bridges else None
    used = [t for t in TYPE_STYLE if any(scene.contacts[k]['type'] == t for k in shown)]
    data = dict(
        title=title or ens.title, kind=ens.kind, chains=[a, b], sides=list(ens.sides), xlabel=ens.xlabel, x=r2(ens.x),
        unit=ens.xlabel.split('(')[-1].rstrip(')'), nframes=ens.n_frames, natoms=len(atoms),
        npoints=scene.points.shape[1],
        atoms=dict(chain=atoms.chain.tolist(), resn=atoms.resname.tolist(), resi=atoms.resnum.tolist(),
                   name=atoms.name.tolist(), elem=atoms.element.tolist()),
        coords=base64.b64encode(coords.tobytes()).decode(),
        types=[dict(name=t, color=_style(t)[0], dash=DASH[_style(t)[1]], badge=_style(t)[2],
                    directional=t in DIRECTIONAL) for t in used],
        contacts=[dict(type=scene.contacts[k]['type'], ra=scene.contacts[k]['residue_a'],
                       rb=scene.contacts[k]['residue_b'], drawn=k in drawn, occ=r2(scene.occ[k]),
                       dist=r2(scene.dist[k]), ia=scene.atom_a[k].tolist(), ib=scene.atom_b[k].tolist(),
                       water=water(k))
                  for k in shown],
        left=[dict(res=r, fill=fill(r), area=r2(scene.area_of(a, r))) for r in scene.left],
        right=[dict(res=r, fill=fill(r), area=r2(scene.area_of(b, r))) for r in scene.right],
        interface={a: sorted({scene.number[(a, scene.contacts[k]['residue_a'])] for k in shown}),
                   b: sorted({scene.number[(b, scene.contacts[k]['residue_b'])] for k in shown})},
        number={c: {r: scene.number[(c, r)] for r in {scene.contacts[k][side] for k in shown}}
                for c, side in ((a, 'residue_a'), (b, 'residue_b'))},
        amino=sorted(r for r in RESIDUE_CLASS if len(r) == 3),
        total=r2(scene.total), bsa=r2(scene.bsa), buriedMin=BURIED_MIN,
        hbonds=r2(scene.type_counts.get('Hydrogen bond', np.zeros(ens.n_frames))),
        curve=r2(ens.curve) if ens.curve is not None else None, curveLabel=ens.curve_label,
    )
    pkg = resources.files('ppitango')
    html = (pkg / 'viewer.html').read_text(encoding='utf-8')
    lib = (pkg / 'vendor' / '3Dmol-min.js').read_text(encoding='utf-8')
    payload = json.dumps(data, separators=(',', ':')).replace('</', '<\\/')
    # str.replace would re-scan inserted text; split keeps the three parts independent
    head, rest = html.split('__3DMOL__', 1)
    mid, tail = rest.split('__DATA__', 1)
    try:
        with open(output, 'w', encoding='utf-8') as f:
            f.write(head.replace('__TITLE__', data['title']) + lib.replace('</script', '<\\/script') + mid + payload
                    + tail)
    except OSError as e:     # a full disk fails in the write, where the error no longer names the file
        raise OSError(e.errno, e.strerror, output) from None
    return output
