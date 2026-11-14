"""Run PandaPPI on every frame and turn the result into smooth per-contact time series."""
import contextlib
import io
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .ensembles import Ensemble

SKIP = ('Van der Waals',)   # atom-level contacts are counted, not drawn
BASIC = {'ARG': ('NH1', 'NH2', 'NE'), 'LYS': ('NZ',), 'HIS': ('ND1', 'NE2')}
PHOSPHATE = ('OP1', 'OP2')


def _phosphate_bridges(pdb_text: str, a: str, b: str, cutoff: float = 4.0):
    """Salt bridges between basic side chains and backbone phosphates (PandaPPI pairs amino acids only)."""
    basic, phos = [], []
    for line in pdb_text.splitlines():
        if not line.startswith('ATOM'):
            continue
        name, resname, chain = line[12:16].strip(), line[17:20].strip(), line[21]
        item = (chain, f'{resname} {int(line[22:26])}', name, (float(line[30:38]), float(line[38:46]), float(line[46:54])))
        if name in BASIC.get(resname, ()):
            basic.append(item)
        elif name in PHOSPHATE:
            phos.append(item)
    best = {}
    for c1, r1, n1, p1 in basic:
        for c2, r2, n2, p2 in phos:
            if c1 == c2:
                continue
            d = sum((u - v) ** 2 for u, v in zip(p1, p2)) ** 0.5
            if d <= cutoff:
                row = (r1, n1, r2, n2) if c1 == a else (r2, n2, r1, n1)
                key = (row[0], row[2])
                if key not in best or d < best[key][-1]:
                    best[key] = ('Salt bridge / ionic', *row, round(d, 2))
    return list(best.values())


def _frame(job):
    """Contacts and buried surface of one frame (runs in a worker process)."""
    pdb_text, a, b, cutoff = job
    from pandappi import PandaPPI
    from pandappi.reports import generator
    with tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False) as f:
        f.write(pdb_text)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            ppi = PandaPPI(f.name, [a, b])
            ppi.map_interactions(distance_cutoff=cutoff)
            df = generator.create_interaction_report(ppi.interactions)
            buried = ppi.buried_surface(a, b)
    finally:
        os.unlink(f.name)
    rows = []
    if df is not None:
        for r in df.itertuples(index=False):
            t, c1, r1, a1, c2, r2, a2, d = r
            if c1 == c2:
                continue
            if c1 != a:  # chain a always on the left
                r1, a1, r2, a2 = r2, a2, r1, a1
            rows.append((t, r1, a1, r2, a2, d))
    rows += _phosphate_bridges(pdb_text, a, b)
    return rows, [(c, r, area) for (c, r), area in buried.items()]


def _smooth(m: np.ndarray, window: int, periodic: bool) -> np.ndarray:
    """Centred rolling mean along the last axis, ignoring NaN."""
    if window <= 1:
        return m
    half = window // 2
    ok = ~np.isnan(m)
    vals = np.where(ok, m, 0.0)
    if periodic:
        idx = (np.arange(m.shape[-1])[:, None] + np.arange(-half, half + 1)[None]) % m.shape[-1]
        total, count = vals[..., idx].sum(-1), ok[..., idx].sum(-1)
    else:
        csum = np.concatenate([np.zeros(m.shape[:-1] + (1,)), np.cumsum(vals, -1)], -1)
        ccnt = np.concatenate([np.zeros(m.shape[:-1] + (1,)), np.cumsum(ok, -1)], -1)
        lo = np.clip(np.arange(m.shape[-1]) - half, 0, None)
        hi = np.clip(np.arange(m.shape[-1]) + half + 1, None, m.shape[-1])
        total, count = csum[..., hi] - csum[..., lo], ccnt[..., hi] - ccnt[..., lo]
    return np.where(count > 0, total / np.maximum(count, 1), np.nan)


@dataclass
class Analysis:
    ensemble: Ensemble
    contacts: pd.DataFrame   # frame, x, type, residue_a, atom_a, residue_b, atom_b, distance
    buried: pd.DataFrame     # frame, x, chain, residue, area
    vdw: np.ndarray          # van der Waals atom contacts per frame
    cutoff: float

    def save(self, prefix: str) -> List[str]:
        files = [f'{prefix}_contacts.csv', f'{prefix}_buried_surface.csv']
        self.contacts.to_csv(files[0], index=False)
        self.buried.to_csv(files[1], index=False)
        return files

    def scene(self, smooth: Optional[int] = None, max_residues: int = 14, min_occupancy: float = 0.5) -> 'Scene':
        return Scene(self, smooth, max_residues, min_occupancy)


def analyze(ensemble: Ensemble, cutoff: float = 4.5, workers: Optional[int] = None,
            progress: bool = True) -> Analysis:
    """PandaPPI contacts and buried surface for every frame, in parallel."""
    a, b = ensemble.chains
    jobs = [(ensemble.pdb(i), a, b, cutoff) for i in range(ensemble.n_frames)]
    workers = workers or min(os.cpu_count() or 1, len(jobs))
    results = []
    with ProcessPoolExecutor(workers) as pool:
        for i, res in enumerate(pool.map(_frame, jobs), 1):
            results.append(res)
            if progress:
                print(f'\r  analysing frames with PandaPPI: {i}/{len(jobs)}', end='', flush=True)
    if progress:
        print()
    contacts, buried, vdw = [], [], np.zeros(len(jobs), int)
    for i, (rows, areas) in enumerate(results):
        x = ensemble.x[i]
        vdw[i] = sum(r[0] == 'Van der Waals' for r in rows)
        contacts += [(i, x, *r) for r in rows if r[0] not in SKIP]
        buried += [(i, x, *r) for r in areas]
    return Analysis(
        ensemble,
        pd.DataFrame(contacts, columns=['frame', 'x', 'type', 'residue_a', 'atom_a', 'residue_b', 'atom_b',
                                        'distance']),
        pd.DataFrame(buried, columns=['frame', 'x', 'chain', 'residue', 'area']), vdw, cutoff)


class Scene:
    """Everything the renderers draw, as arrays over frames.

    contacts   one entry per (type, residue_a, residue_b) ever seen
    occ        (contacts, frames) smoothed presence, 0..1: the opacity of a contact line
    dist       (contacts, frames) distance between the contact's two atoms, measured in every frame
    atom_a/b   (contacts, frames) indices into ensemble.atoms of those two atoms
    drawn      indices of the contacts shown on the 2D map; left/right its residues, top to bottom
    """

    def __init__(self, analysis: Analysis, smooth: Optional[int], max_residues: int, min_occupancy: float):
        ens = analysis.ensemble
        self.analysis, self.ensemble = analysis, ens
        F = ens.n_frames
        if smooth is None:
            smooth = 1 if ens.kind == 'nma' else 5
        self.smooth = smooth
        a, b = ens.chains
        atoms = ens.atoms
        label = atoms.resname + ' ' + atoms.resnum.astype(str)
        by_name = {k: i for i, k in enumerate(zip(atoms.chain, label, atoms.name))}
        by_res: Dict[tuple, np.ndarray] = {k: g.values for k, g in
                                           pd.Series(atoms.index).groupby([atoms.chain, label])}

        df = analysis.contacts
        # the closest atom pair of each contact in each frame
        df = df.sort_values('distance').drop_duplicates(['frame', 'type', 'residue_a', 'residue_b'])
        keys = list(df.groupby(['type', 'residue_a', 'residue_b'], sort=False).size().index)
        kid = {k: i for i, k in enumerate(keys)}
        K = len(keys)
        present = np.zeros((K, F))
        ea, eb = -np.ones((K, F), int), -np.ones((K, F), int)
        for r in df.itertuples(index=False):
            k = kid[(r.type, r.residue_a, r.residue_b)]
            present[k, r.frame] = 1
            i, j = by_name.get((a, r.residue_a, r.atom_a)), by_name.get((b, r.residue_b, r.atom_b))
            if i is None or j is None:  # ring centroids and the like: use the closest atoms instead
                ia, ib = by_res[(a, r.residue_a)], by_res[(b, r.residue_b)]
                d = np.linalg.norm(ens.coords[r.frame, ia][:, None] - ens.coords[r.frame, ib][None], axis=-1)
                p, q = np.unravel_index(d.argmin(), d.shape)
                i, j = ia[p], ib[q]
            ea[k, r.frame], eb[k, r.frame] = i, j
        # while a contact is absent keep pointing at the atoms it last used, so its line fades in place
        frames = np.arange(F)
        for k in range(K):
            seen = frames[ea[k] >= 0]
            nearest = seen[np.abs(frames[:, None] - seen[None]).argmin(1)]
            ea[k], eb[k] = ea[k, nearest], eb[k, nearest]
        self.contacts = [dict(type=t, residue_a=ra, residue_b=rb) for t, ra, rb in keys]
        self.atom_a, self.atom_b = ea, eb
        self.dist = np.linalg.norm(ens.coords[frames[None], ea] - ens.coords[frames[None], eb], axis=-1) \
            if K else np.zeros((0, F))
        self.present = present
        self.occ = _smooth(present, smooth, ens.periodic)

        types = [c['type'] for c in self.contacts]
        self.types = sorted(set(types), key=lambda t: -self.occ[[x == t for x in types]].sum())
        self.type_counts = {t: self.occ[[x == t for x in types]].sum(0) for t in self.types}
        self.total = self.occ.sum(0) if K else np.zeros(F)
        self.vdw = _smooth(analysis.vdw[None].astype(float), smooth, ens.periodic)[0]

        # buried surface per residue
        bur = analysis.buried
        self.residues = list(bur.groupby(['chain', 'residue'], sort=False).size().index)
        rid = {k: i for i, k in enumerate(self.residues)}
        area = np.zeros((len(self.residues), F))
        for r in bur.itertuples(index=False):
            area[rid[(r.chain, r.residue)], r.frame] = r.area
        self.area = _smooth(area, smooth, ens.periodic)
        self.bsa = self.area.sum(0)
        self.bsa_chain = {c: self.area[[k[0] == c for k in self.residues]].sum(0) for c in (a, b)}

        # what goes on the 2D map: contacts that are really there at some point, on the busiest residues
        weight = self.occ.mean(1)
        cand = [k for k in range(K) if self.occ[k].max() >= min_occupancy]
        left = self._busiest([self.contacts[k]['residue_a'] for k in cand], weight[cand], max_residues)
        right = self._busiest([self.contacts[k]['residue_b'] for k in cand], weight[cand], max_residues)
        self.drawn = [k for k in cand if self.contacts[k]['residue_a'] in left
                      and self.contacts[k]['residue_b'] in right]
        num = lambda res: int(res.split()[1])
        self.left = sorted({self.contacts[k]['residue_a'] for k in self.drawn}, key=num)
        lpos = {r: i for i, r in enumerate(self.left)}
        partners: Dict[str, list] = {}
        for k in self.drawn:
            c = self.contacts[k]
            partners.setdefault(c['residue_b'], []).append((lpos[c['residue_a']], weight[k]))
        # order the right side by where its partners sit, to cut line crossings
        self.right = sorted(partners, key=lambda r: (np.average([p for p, _ in partners[r]],
                                                                weights=[w for _, w in partners[r]]), num(r)))
        self.hidden = len(cand) - len(self.drawn)

    @staticmethod
    def _busiest(residues, weights, n) -> set:
        score: Dict[str, float] = {}
        for r, w in zip(residues, weights):
            score[r] = score.get(r, 0) + w
        return set(sorted(score, key=score.get, reverse=True)[:n])

    def area_of(self, chain: str, residue: str) -> np.ndarray:
        try:
            return self.area[self.residues.index((chain, residue))]
        except ValueError:
            return np.zeros(self.ensemble.n_frames)
