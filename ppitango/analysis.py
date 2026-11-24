"""Run PandaPPI on every frame and turn the result into smooth per-contact time series."""
import contextlib
import io
import os
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .ensembles import Ensemble

SKIP = ('Van der Waals',)   # atom-level contacts are counted, not drawn
BASIC = {'ARG': ('NH1', 'NH2', 'NE'), 'LYS': ('NZ',), 'HIS': ('ND1', 'NE2')}
CATION = {'ARG': 'CZ', 'LYS': 'NZ'}
SIX, FIVE = ('N1', 'C2', 'N3', 'C4', 'C5', 'C6'), ('C4', 'C5', 'N7', 'C8', 'N9')
RINGS = {'PHE': [('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')], 'TYR': [('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')],
         'TRP': [('CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'), ('CG', 'CD1', 'NE1', 'CE2', 'CD2')],
         'HIS': [('CG', 'ND1', 'CD2', 'CE1', 'NE2')],
         **{r: [SIX, FIVE] for r in ('DA', 'DG', 'A', 'G')}, **{r: [SIX] for r in ('DC', 'DT', 'DU', 'C', 'U')}}
AMIDES = {'ASN': [('ND2', 'CG', 'OD1')], 'GLN': [('NE2', 'CD', 'OE1')]}   # besides the backbone amide (N, C, O)
BASES = ('DA', 'DG', 'DC', 'DT', 'DU', 'A', 'G', 'C', 'U')
AMINO = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
         'SER', 'THR', 'TRP', 'TYR', 'VAL'}
GREASY = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR', 'CYS'}
METALS = {'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'ZN', 'MN', 'FE', 'CU', 'NI', 'CO', 'CD', 'SR', 'BA', 'HG'}
HALIDES = {'F', 'CL', 'BR', 'I'}


def _extra_contacts(pdb_text: str, a: str, b: str):
    """Contacts PandaPPI does not look for, since it pairs amino acids with amino acids. With its cutoffs:

    - salt bridges from Arg, Lys, His to backbone phosphates (4.0 Å)
    - stacking of an aromatic side chain on a base: ring centres within 7.0 Å, planes parallel or at a right
      angle (within 30°)
    - cation-π from Arg or Lys to a base (6.0 Å to the ring centre)
    - for a ligand, sugar, lipid or ion on one side: hydrophobic contacts between its carbons that carry no
      N or O and side-chain carbons of hydrophobic residues (5.0 Å), metal ions on N, O, S (3.0 Å), and halide
      ions on basic groups (4.0 Å)
    """
    rows = []
    for line in pdb_text.splitlines():
        if line.startswith('ATOM'):
            rows.append((line[21], line[17:20].strip(), int(line[22:26]), line[12:16].strip(),
                         line[76:78].strip().upper(), float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not rows:
        return []
    chain, resname, resnum, name, element = (np.array(c) for c in list(zip(*rows))[:5])
    xyz = np.array([r[5:] for r in rows])
    res = np.array([f'{r} {n}' for r, n in zip(resname, resnum)])
    polymer = np.isin(resname, list(AMINO)) | np.isin(resname, BASES)
    side = {a: chain == a, b: chain == b}
    best = {}

    def add(kind, i, j, d, centre=()):   # i, j: atom rows on either side; keep the closest pair per residue pair
        # `centre`: the rows that stand for a ring. The distance is to its centre, so no atom is named there
        named = ['' if k in centre else name[k] for k in (i, j)]
        if chain[i] != a:
            i, j, named = j, i, named[::-1]
        key = (kind, res[i], res[j])
        if key not in best or d < best[key][-1]:
            best[key] = (kind, res[i], named[0], res[j], named[1], round(float(d), 2))

    def pairs(kind, first, second, cutoff):
        for one, two in ((a, b), (b, a)):
            i, j = np.where(first & side[one])[0], np.where(second & side[two])[0]
            if len(i) and len(j):
                d = np.linalg.norm(xyz[i][:, None] - xyz[j][None], axis=-1)
                for p, q in zip(*np.where(d <= cutoff)):
                    add(kind, i[p], j[q], d[p, q])

    basic = np.array([n in BASIC.get(r, ()) for r, n in zip(resname, name)])
    pairs('Salt bridge / ionic', basic, np.isin(name, ('OP1', 'OP2')), 4.0)

    rings = []                       # (chain, is base, first atom row, centre, normal)
    start = np.r_[0, np.where((chain[1:] != chain[:-1]) | (res[1:] != res[:-1]))[0] + 1, len(res)]
    for s, e in zip(start[:-1], start[1:]):
        for ring in RINGS.get(resname[s], ()):
            at = [s + k for k in range(e - s) if name[s + k] in ring]
            if len(at) == len(ring):
                centre = xyz[at].mean(0)
                rings.append((chain[s], resname[s] in BASES, at[0], centre, np.linalg.svd(xyz[at] - centre)[2][2]))
    for c1, base1, i, p1, n1 in rings:
        for c2, base2, j, p2, n2 in rings:
            if c1 == a and c2 == b and (base1 != base2):
                d = np.linalg.norm(p1 - p2)
                angle = np.degrees(np.arccos(min(abs(float(n1 @ n2)), 1.0)))
                if d <= 7.0 and (angle <= 30 or angle >= 60):
                    add('Pi-pi stacking', i, j, d, centre=(i, j))
    for is_cation in (np.array([CATION.get(r) == n for r, n in zip(resname, name)]),):
        for c, base, j, centre, _ in rings:
            if base:
                i = np.where(is_cation & ~side[c] & (side[a] | side[b]))[0]
                for k in i[np.linalg.norm(xyz[i] - centre, axis=1) <= 6.0] if len(i) else ():
                    add('Cation-pi', k, j, np.linalg.norm(xyz[k] - centre), centre=(j,))

    if not polymer.all():
        other = ~polymer
        lone = np.array([(res == r).sum() == 1 for r in res]) & other
        metal, halide = lone & np.isin(element, list(METALS)), lone & np.isin(element, list(HALIDES))
        pairs('Metal coordination', metal, np.isin(element, ('N', 'O', 'S')) & ~lone, 3.0)
        pairs('Salt bridge / ionic', halide, basic, 4.0)
        carbon = other & (element == 'C')
        hetero = np.where(other & np.isin(element, ('N', 'O')))[0]
        if carbon.any() and len(hetero):
            idx = np.where(carbon)[0]
            polar = (np.linalg.norm(xyz[idx][:, None] - xyz[hetero][None], axis=-1) < 1.7).any(1)
            carbon[idx[polar]] = False
        greasy = np.isin(resname, list(GREASY)) & (element == 'C') & ~np.isin(name, ('C', 'CA'))
        pairs('Hydrophobic', carbon, greasy | carbon, 5.0)
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
    rows += _extra_contacts(pdb_text, a, b)
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
        for table, file in zip((self.contacts, self.buried), files):
            try:
                table.to_csv(file, index=False)
            except OSError as e:     # a full disk fails in the write, where the error no longer names the file
                raise OSError(e.errno, e.strerror, file) from None
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
            if progress and sys.stdout.isatty():   # a counter in a terminal, one line in a log file
                print(f'\r  analysing frames with PandaPPI: {i}/{len(jobs)}', end='', flush=True)
    if progress:
        print('' if sys.stdout.isatty() else f'  analysed {len(jobs)} frames with PandaPPI')
    # PandaPPI names a residue by what the PDB frame holds (three letters, the frame's own number): show it by
    # its full name and original number instead
    at = ensemble.atoms
    shown = {(c, f'{r[:3]} {n}'): label for c, r, n, label in zip(at.chain, at.resname, at.resnum, ensemble.labels)}
    contacts, buried, vdw = [], [], np.zeros(len(jobs), int)
    for i, (rows, areas) in enumerate(results):
        x = ensemble.x[i]
        vdw[i] = sum(r[0] == 'Van der Waals' for r in rows)
        contacts += [(i, x, t, shown.get((a, r1), r1), a1, shown.get((b, r2), r2), a2, d)
                     for t, r1, a1, r2, a2, d in rows if t not in SKIP]
        buried += [(i, x, c, shown.get((c, r), r), area) for c, r, area in areas]
    return Analysis(
        ensemble,
        pd.DataFrame(contacts, columns=['frame', 'x', 'type', 'residue_a', 'atom_a', 'residue_b', 'atom_b',
                                        'distance']),
        pd.DataFrame(buried, columns=['frame', 'x', 'chain', 'residue', 'area']), vdw, cutoff)


class Scene:
    """Everything the renderers draw, as arrays over frames.

    contacts   one entry per (type, residue_a, residue_b) ever seen
    occ        (contacts, frames) smoothed presence, 0..1: the opacity of a contact line
    points     (frames, atoms + centres, 3) the atoms of the ensemble, then the centres of the rings and amide
               groups that contacts are measured to (cation-π, CH-π, π-π stacking, amide-aromatic)
    dist       (contacts, frames) distance between the contact's two end points, measured in every frame
    atom_a/b   (contacts, frames) indices into `points` of those two ends: an atom, or the centre of a group
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
        label = ens.labels
        self.rank = {k: i for i, k in reversed(list(enumerate(zip(atoms.chain, label))))}     # order along a chain
        self.number = dict(zip(zip(atoms.chain, label), atoms.resnum.tolist()))   # the number in the PDB frames
        by_name = {k: i for i, k in enumerate(zip(atoms.chain, label, atoms.name))}
        by_res: Dict[tuple, np.ndarray] = {k: g.values for k, g in
                                           pd.Series(atoms.index).groupby([atoms.chain, label])}
        N = len(atoms)
        names, resnames = atoms.name.values, atoms.resname.values
        groups: Dict[tuple, int] = {}      # atoms of a ring or amide group -> index of its centre in `points`

        def ends(chain, residue, atom):
            """What a contact may end on at this residue: its named atom, or else its rings and amide groups."""
            i = by_name.get((chain, residue, atom))
            if i is not None:
                return [(i,)]
            rows = by_res[(chain, residue)]
            found = []
            for group in RINGS.get(resnames[rows[0]], []) + AMIDES.get(resnames[rows[0]], []) + [('N', 'C', 'O')]:
                at = tuple(rows[np.isin(names[rows], group)])
                if len(at) == len(group):
                    found.append(at)
            return found

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
            # a contact to a ring or an amide group names no atom on that side: its distance is to the centre of
            # the group. Which group it was (a tryptophan has two rings) is the one that gives that distance
            xyz = ens.coords[r.frame]
            best = None
            for ga in ends(a, r.residue_a, r.atom_a):
                for gb in ends(b, r.residue_b, r.atom_b):
                    off = abs(np.linalg.norm(xyz[list(ga)].mean(0) - xyz[list(gb)].mean(0)) - r.distance)
                    if best is None or off < best[0]:
                        best = (off, ga, gb)
            if best is not None and best[0] < 0.05:
                i, j = (g[0] if len(g) == 1 else groups.setdefault(g, N + len(groups)) for g in best[1:])
            else:                       # nothing known fits: the closest atoms of the two residues
                ia, ib = by_res[(a, r.residue_a)], by_res[(b, r.residue_b)]
                d = np.linalg.norm(xyz[ia][:, None] - xyz[ib][None], axis=-1)
                p, q = np.unravel_index(d.argmin(), d.shape)
                i, j = ia[p], ib[q]
            ea[k, r.frame], eb[k, r.frame] = i, j
        self.points = np.concatenate([ens.coords] + [ens.coords[:, list(g)].mean(1, keepdims=True) for g in groups],
                                     axis=1)
        self.point_chain = np.r_[atoms.chain.values, [atoms.chain.values[g[0]] for g in groups]]
        # while a contact is absent keep pointing at the atoms it last used, so its line fades in place
        frames = np.arange(F)
        for k in range(K):
            seen = frames[ea[k] >= 0]
            nearest = seen[np.abs(frames[:, None] - seen[None]).argmin(1)]
            ea[k], eb[k] = ea[k, nearest], eb[k, nearest]
        self.contacts = [dict(type=t, residue_a=ra, residue_b=rb) for t, ra, rb in keys]
        self.atom_a, self.atom_b = ea, eb
        self.dist = np.linalg.norm(self.points[frames[None], ea] - self.points[frames[None], eb], axis=-1) \
            if K else np.zeros((0, F))
        self.present = present
        # the water each water-mediated contact goes through, per frame: the one closest to both of its atoms
        self.bridge = {}
        if ens.waters is not None:
            for k, (t, _, _) in enumerate(keys):
                if t != 'Water-mediated':
                    continue
                pos = np.full((F, 3), np.nan)
                for f in np.where(present[k] > 0)[0]:
                    w = ens.waters[f]
                    if len(w):
                        reach = np.maximum(np.linalg.norm(w - self.points[f, ea[k, f]], axis=1),
                                           np.linalg.norm(w - self.points[f, eb[k, f]], axis=1))
                        if reach.min() <= 4.5:
                            pos[f] = w[reach.argmin()]
                self.bridge[k] = pos
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
        self.left = sorted({self.contacts[k]['residue_a'] for k in self.drawn}, key=lambda r: self.rank[(a, r)])
        lpos = {r: i for i, r in enumerate(self.left)}
        partners: Dict[str, list] = {}
        for k in self.drawn:
            c = self.contacts[k]
            partners.setdefault(c['residue_b'], []).append((lpos[c['residue_a']], weight[k]))
        # order the right side by where its partners sit, to cut line crossings
        self.right = sorted(partners, key=lambda r: (np.average([p for p, _ in partners[r]],
                                                                weights=[w for _, w in partners[r]]), self.rank[(b, r)]))
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
