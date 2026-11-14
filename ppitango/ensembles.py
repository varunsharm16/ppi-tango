"""Where the frames come from: a normal mode, an MD trajectory, or FEP lambda windows.

Every source returns an Ensemble: the heavy atoms of two chains, their coordinates in each frame,
and the axis the frames run along (mode phase, time, lambda).
"""
from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

# protonation / force-field variants -> the names PandaPPI's detectors know
RESNAME_ALIASES = {'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS', 'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
                   'CYX': 'CYS', 'CYM': 'CYS', 'ASH': 'ASP', 'GLH': 'GLU', 'LYN': 'LYS', 'MSE': 'MET'}
# force-field atom names -> PDB names (terminal oxygens, CHARMM's isoleucine CD)
ATOM_ALIASES = {('ILE', 'CD'): 'CD1', 'OT1': 'O', 'OC1': 'O', 'O1': 'O', 'OT2': 'OXT', 'OC2': 'OXT', 'O2': 'OXT'}
# nucleotides: CHARMM words, Amber terminal forms (DA5, DG3, DCN) and old PDB spellings -> DA DT DG DC / A U G C
NUCLEOTIDES = {'ADE': 'A', 'GUA': 'G', 'CYT': 'C', 'THY': 'T', 'URA': 'U'}
WATER = {'HOH', 'WAT', 'SOL', 'TIP3', 'TIP', 'TIP4', 'TP3', 'T3P', 'T4P', 'SPC', 'OPC', 'H2O'}
SOLVENT = {'HOH', 'WAT', 'SOL', 'TIP3', 'TIP', 'NA', 'CL', 'K', 'NA+', 'CL-', 'SOD', 'CLA', 'POT', 'MG', 'CA', 'ZN'}


@dataclass
class Ensemble:
    atoms: pd.DataFrame            # chain, resname, resnum, name, element: one row per heavy atom
    coords: np.ndarray             # (frames, atoms, 3) in Å
    chains: tuple                  # (chain_a, chain_b)
    x: np.ndarray                  # axis value of each frame
    xlabel: str
    kind: str                      # 'nma' | 'md' | 'fep'
    title: str = ''
    periodic: bool = False         # last frame flows back into the first (a mode oscillation)
    curve: Optional[np.ndarray] = None   # an extra per-frame series drawn on the timeline (FEP: ΔG)
    curve_label: str = ''
    notes: List[str] = field(default_factory=list)
    waters: Optional[list] = None  # per frame: oxygen positions of the waters between the two sides
    log: List[str] = field(default_factory=list)     # what was decided while reading (chains found, box handling)

    @property
    def n_frames(self) -> int:
        return len(self.coords)

    def pdb(self, frame: int) -> str:
        """One frame as PDB text."""
        return _pdb_text(self.atoms, self.coords[frame], self.waters[frame] if self.waters is not None else None)


def _pdb_text(atoms: pd.DataFrame, xyz: np.ndarray, waters: Optional[np.ndarray] = None) -> str:
    lines = []
    for i, (a, (x, y, z)) in enumerate(zip(atoms.itertuples(), xyz), 1):
        name = a.name if len(a.name) == 4 else f' {a.name:<3}'
        lines.append(f'ATOM  {i % 100000:5d} {name} {a.resname:>3} {a.chain}{a.resnum:4d}    '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.element:>2}')
    if waters is not None:   # interface waters ride along in the first chain, as PDB files carry them
        chain = atoms.chain.iloc[0]
        for k, (x, y, z) in enumerate(waters, 1):
            lines.append(f'HETATM{(len(atoms) + k) % 100000:5d}  O   HOH {chain}{k % 10000:4d}    '
                         f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O')
    return '\n'.join(lines) + '\nEND\n'


def canonical(resname: str, name: str, ribose: bool = False):
    """Force-field residue and atom names -> the PDB names PandaPPI's detectors know."""
    name = name.replace('*', "'")
    name = {'O1P': 'OP1', 'O2P': 'OP2', 'C5M': 'C7'}.get(name, name)
    resname = RESNAME_ALIASES.get(resname, resname)
    base = resname
    if len(base) == 3 and base[0] in 'DR' and base[1] in 'ACGTU' and base[2] in '53N':   # Amber DA5, RU3, DCN
        base = base[:2]
    elif len(base) == 2 and base[0] in 'ACGTU' and base[1] in '53':                      # Amber RNA A5, U3
        base = base[0]
    if base in NUCLEOTIDES:                         # CHARMM uses one word for DNA and RNA: the sugar decides
        base = NUCLEOTIDES[base] if ribose else 'D' + NUCLEOTIDES[base]
    elif base in ('RA', 'RC', 'RG', 'RU'):
        base = base[1]
    if base in ('DA', 'DT', 'DG', 'DC', 'DU', 'A', 'U', 'G', 'C'):
        return base, name
    name = ATOM_ALIASES.get((resname, name), ATOM_ALIASES.get(name, name))
    return resname, name


def read_pdb(path: str, chains: Optional[Sequence[str]] = None, all_models: bool = False):
    """Heavy polymer atoms of a PDB file.

    Returns (atoms, coords, index): coords is (models, atoms, 3) and index gives each kept atom's
    position among all ATOM/HETATM records of a model, which is how trajectory frames are sliced.
    """
    rows, index, models, xyz, n = [], [], [], [], -1
    first = True
    with open(path) as f:
        for line in f:
            rec = line[:6]
            if rec.startswith('ENDMDL'):
                models.append(xyz)
                xyz, n, first = [], -1, False
                if not all_models:
                    break
                continue
            if rec not in ('ATOM  ', 'HETATM'):
                continue
            n += 1
            name, resname, chain = line[12:16].strip(), line[17:20].strip(), line[21]
            element = (line[76:78].strip() or name.lstrip('0123456789')[:1]).upper()
            if (rec == 'HETATM' and resname != 'MSE') or resname in SOLVENT or element in ('H', 'D') \
                    or line[16] not in ' A' or (chains and chain not in chains):
                continue
            if first:
                resname, name = canonical(resname, name)
                rows.append((chain, resname, int(line[22:26]), name, element))
                index.append(n)
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if xyz:
        models.append(xyz)
    atoms = pd.DataFrame(rows, columns=['chain', 'resname', 'resnum', 'name', 'element'])
    if atoms.empty:
        raise ValueError(f'No polymer atoms found in {path}' + (f' for chains {list(chains)}' if chains else ''))
    if len({len(m) for m in models}) != 1:
        raise ValueError(f'{path}: models do not all contain the same atoms')
    return atoms, np.array(models, dtype=float), np.array(index)


def pick_chains(path: str, chains: Optional[Sequence[str]]) -> tuple:
    """The two chains to follow: the ones asked for, or the pair with the most atoms in contact."""
    if chains:
        if len(chains) != 2:
            raise ValueError('Give exactly two chains, e.g. --chains A D')
        return tuple(chains)
    atoms, coords, _ = read_pdb(path)
    ids = list(dict.fromkeys(atoms.chain))
    if len(ids) < 2:
        raise ValueError(f'{path} has a single chain ({ids}); an interface needs two')
    best, pair = -1, None
    xyz = {c: coords[0][(atoms.chain == c).values] for c in ids}
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            # atoms of b within 4.5 Å of a, one atom of b at a time to keep memory flat
            n = sum(((xyz[a] - p) ** 2).sum(1).min() < 4.5 ** 2 for p in xyz[b])
            if n > best:
                best, pair = n, (a, b)
    return pair


def _superpose(coords: np.ndarray, sel: np.ndarray, extras: Optional[list] = None) -> np.ndarray:
    """Fit every frame onto the first, using the selected atoms (Kabsch). `extras` holds further points per
    frame (waters) that are moved along, in place."""
    ref = coords[0, sel]
    ref_c = ref.mean(0)
    out = np.empty_like(coords)
    for i, xyz in enumerate(coords):
        c = xyz[sel].mean(0)
        u, _, vt = np.linalg.svd((xyz[sel] - c).T @ (ref - ref_c))
        d = np.sign(np.linalg.det(u @ vt))
        rot = u @ np.diag([1, 1, d]) @ vt
        out[i] = (xyz - c) @ rot + ref_c
        if extras is not None:
            extras[i] = (extras[i] - c) @ rot + ref_c
    return out


def _name(path: str) -> str:
    import os
    return os.path.splitext(os.path.basename(path))[0]


# ---------------------------------------------------------------- normal modes

def anm_modes(ca: np.ndarray, cutoff: float = 15.0, n_modes: int = 20):
    """Anisotropic network model on Cα atoms. Returns (eigenvalues, modes[n_modes, N, 3]), slowest first."""
    n = len(ca)
    d = ca[:, None, :] - ca[None, :, :]
    r2 = (d ** 2).sum(-1)
    k = np.where((r2 < cutoff ** 2) & (r2 > 0), -1.0 / np.where(r2 > 0, r2, 1), 0.0)
    blocks = k[:, :, None, None] * d[:, :, :, None] * d[:, :, None, :]   # off-diagonal 3x3 blocks
    hess = blocks.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    diag = -blocks.sum(1)
    for i in range(n):
        hess[3 * i:3 * i + 3, 3 * i:3 * i + 3] = diag[i]
    vals, vecs = np.linalg.eigh(hess)
    if vals[6] < 1e-6:
        raise ValueError('The elastic network is disconnected (more than 6 zero modes): '
                         'the two chains do not touch within the cutoff')
    return vals[6:6 + n_modes], vecs[:, 6:6 + n_modes].T.reshape(-1, n, 3)


def nma(pdb: str, chains: Optional[Sequence[str]] = None, mode: int = 1, rmsd: float = 1.0,
        frames: int = 48, cutoff: float = 15.0, select: Optional[Sequence[str]] = None,
        renumber: Optional[dict] = None) -> Ensemble:
    """One oscillation along an elastic-network normal mode of the two-sided complex.

    The network has one node per residue: Cα for amino acids, C4' for nucleotides. Each residue moves rigidly
    with its node, and the mode is scaled so the node RMSD from the starting structure peaks at `rmsd` Å.
    `pdb` is a PDB file; with `select`, or for any other structure file, it is read through MDAnalysis.
    """
    log = []
    if select or renumber or not pdb.lower().endswith(('.pdb', '.ent')):
        u = _universe(pdb, None)
        atoms, index, (a, b), log = _sides(u, chains, select, renumber)
        xyz = u.atoms.positions[index].astype(float)
        box = u.trajectory.ts.dimensions
        if box is not None and box[:3].min() > 0:
            xyz = _whole(xyz, (atoms.chain == b).values, u.trajectory.ts.triclinic_dimensions.astype(float))
    else:
        a, b = pick_chains(pdb, chains)
        atoms, coords, _ = read_pdb(pdb, (a, b))
        xyz = coords[0]
    res_key = atoms.chain + atoms.resnum.astype(str)
    is_node = atoms.name.isin(['CA', "C4'"]).values
    node_of = dict(zip(res_key[is_node], range(is_node.sum())))
    missing = ~res_key.isin(node_of)
    if missing.any():  # residues without a node atom (truncated residues) cannot be moved
        atoms, xyz, res_key, is_node = atoms[~missing.values].reset_index(drop=True), xyz[~missing.values], \
            res_key[~missing.values], is_node[~missing.values]
    gone = [c for c in (a, b) if not (atoms.chain == c).any()]
    if gone:
        raise ValueError(f"Chain {' and '.join(gone)} has no Cα or C4' atoms to build the elastic network on")
    vals, modes = anm_modes(xyz[is_node], cutoff, n_modes=max(mode, 1))
    disp = modes[mode - 1]
    disp = disp / np.sqrt((disp ** 2).sum(-1).mean()) * rmsd      # node RMSD = rmsd at full swing
    per_atom = disp[res_key.map(node_of).values]
    phase = np.arange(frames) / frames
    stack = xyz[None] + np.sin(2 * np.pi * phase)[:, None, None] * per_atom[None]
    nodes = "Cα" + (" and C4'" if (atoms.name == "C4'").any() else '')
    return Ensemble(atoms, stack, (a, b), x=phase, xlabel=f'Phase of normal mode {mode} (cycles)', kind='nma',
                    title=f'{_name(pdb)} · normal mode {mode}', periodic=True, log=log,
                    notes=[f'Elastic network model ({nodes}, {cutoff:g} Å cutoff), mode {mode}, ±{rmsd:g} Å RMSD'])


# ---------------------------------------------------------------- trajectories

TIMED = ('XTC', 'TRR', 'NCDF', 'TNG', 'H5MD')   # trajectory formats that store a real time (ps) per frame
LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


def _attr(group, name):
    """A topology attribute as an array of strings, or None when the file format does not carry it."""
    try:
        return np.array([str(v).strip() for v in getattr(group, name)])
    except Exception:
        return None


def _universe(top: str, traj: Optional[str]):
    try:
        import MDAnalysis as mda
    except ImportError as e:
        raise ImportError('Reading trajectory files needs MDAnalysis: pip install "ppi-tango[md]"') from e
    import os
    import warnings
    for f in (top, traj):
        if f and not os.path.exists(f):
            raise FileNotFoundError(f'No such file: {f}')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # MDAnalysis is chatty about PDB records and future changes
        try:
            u = mda.Universe(top, traj) if traj else mda.Universe(top)
        except Exception as e:
            hint = ''
            if top.lower().endswith('.top'):
                hint = (' A GROMACS .top is not a structure file: give the .tpr, .gro or a PDB written from the '
                        'same system instead.')
            elif 'tpx version' in str(e):
                hint = (' The installed MDAnalysis is older than the GROMACS that wrote this .tpr: upgrade it '
                        '(pip install -U MDAnalysis) or give the .gro of the same system as topology.')
            elif 'atoms' in str(e).lower() or 'natoms' in str(e).lower():
                hint = (' The topology and the trajectory must hold the same atoms in the same order: if the '
                        'trajectory was saved for a subset (protein only), write the topology for that subset too.')
            raise ValueError(f'Could not read {top}' + (f' with {traj}' if traj else '') +
                             f': {type(e).__name__}: {e}.{hint}') from e
    return u


def _polymer(u, renumber: Optional[dict] = None):
    """Heavy atoms of protein and nucleic acid residues, with a chain label for every atom.

    Returns (atoms, index, table): `atoms` has one row per kept atom (chain, segid, resname, resnum, name,
    element), `index` is each row's position in the trajectory's atom order and `table` describes the chains.
    Residues count as polymer when they have a backbone (N and CA, or a sugar C4'), whatever they are called,
    so force-field names (HIE, CYX, ...) and modified residues stay and water, ions, lipids and ligands go.
    """
    ag = u.atoms
    names = _attr(ag, 'names')
    resix = ag.resindices
    n_res = resix.max() + 1
    has = lambda *wanted: np.bincount(resix[np.isin(names, wanted)], minlength=n_res) > 0
    protein = has('N') & has('CA')
    nucleic = has("C4'", 'C4*')
    elements = _attr(ag, 'elements')
    bare = np.array([n.lstrip('0123456789') for n in names])
    if elements is not None and (elements != '').all():
        hydrogen = np.isin(np.char.upper(elements), ('H', 'D'))
    else:
        hydrogen = np.array([n[:1] == 'H' for n in bare])
        elements = np.array([n[:1] for n in bare])
    altloc = _attr(ag, 'altLocs')
    keep = (protein | nucleic)[resix] & ~hydrogen
    if altloc is not None:
        keep &= np.isin(altloc, ('', 'A'))
    index = np.where(keep)[0]
    if not len(index):
        raise ValueError('No protein or nucleic acid found in the topology')

    ribose = has("O2'", 'O2*')[resix[index]]
    pairs = [canonical(r, n, rib) for r, n, rib in zip(_attr(ag, 'resnames')[index], names[index], ribose)]
    resname, name = np.array([q[0] for q in pairs]), np.array([q[1] for q in pairs])
    resnum = ag.resids[index]
    segid = _attr(ag, 'segids')
    segid = segid[index] if segid is not None else np.full(len(index), '')
    chain_id = _attr(ag, 'chainIDs')
    chain_id = chain_id[index] if chain_id is not None else np.full(len(index), '')

    if len(set(chain_id) - {'', 'X', 'SYSTEM'}) >= 2 and (chain_id != '').all():
        chain, how = chain_id, 'chain IDs of the topology'
    else:
        # no usable chain IDs (GROMACS, Amber, most PSF files): a new chain starts where the segment changes,
        # the residue numbering restarts, or the backbone is not bonded to the previous residue
        xyz = ag.positions[index].astype(float)
        box = u.trajectory.ts.dimensions
        from MDAnalysis.lib.distances import calc_bonds
        res_of = resix[index]
        starts = np.r_[0, np.where(np.diff(res_of) != 0)[0] + 1]          # first atom of every residue
        res_ids = res_of[starts]
        def pick(wanted):   # residue -> row of its atom with one of these names
            rows = np.where(np.isin(name, wanted))[0]
            return dict(zip(res_of[rows], rows))
        tail, head = pick(('C', "O3'", 'O3*')), pick(('N', 'P'))
        label, k = {}, 0
        for prev, cur, s in zip(res_ids[:-1], res_ids[1:], starts[1:]):
            label[prev] = k
            new = segid[s] != segid[s - 1] or resnum[s] < resnum[s - 1] or protein[cur] != protein[prev]
            if not new and prev in tail and cur in head:
                gap = calc_bonds(xyz[tail[prev]][None], xyz[head[cur]][None],
                                 box=box if box is not None and box[:3].min() > 0 else None)[0]
                new = gap > 2.5
            k += bool(new)
        label[res_ids[-1]] = k
        if k + 1 > len(LETTERS):
            raise ValueError(f'{k + 1} separate chains found; choose the two sides with --select')
        chain = np.array([LETTERS[label[r]] for r in res_of])
        how = 'no chain IDs in the topology: chains lettered in file order, split where the backbone breaks'

    atoms = pd.DataFrame({'chain': chain, 'segid': segid, 'resname': resname, 'resnum': resnum, 'name': name,
                          'element': np.char.upper(elements[index])})
    for c, first in (renumber or {}).items():   # tpr and prmtop files count residues 1, 2, 3 ... through all chains
        if c not in set(chain):
            raise ValueError(f'--renumber: no chain {c}. Chains: {", ".join(dict.fromkeys(chain))}')
        rows = atoms.chain == c
        atoms.loc[rows, 'resnum'] += first - atoms.resnum[rows].iloc[0]
    rows = []
    for c, g in atoms.groupby('chain', sort=False):
        first, last = g.iloc[0], g.iloc[-1]
        n = g[['segid', 'resnum']].drop_duplicates().shape[0]
        seg = f', segment {first.segid}' if first.segid and first.segid != 'SYSTEM' else ''
        rows.append(f'  {c}: {n} residues, {first.resname}{first.resnum} to {last.resname}{last.resnum}{seg}')
    return atoms, index, f'Chains ({how}):\n' + '\n'.join(rows)


def _closest_pair(atoms: pd.DataFrame, xyz: np.ndarray, box) -> tuple:
    """The two chains with the most atoms in contact (4.5 Å), counting across the periodic box."""
    from MDAnalysis.lib.distances import capped_distance
    ids = list(dict.fromkeys(atoms.chain))
    if len(ids) < 2:
        raise ValueError(f'Only one chain found ({ids[0]}); an interface needs two')
    box = box if box is not None and box[:3].min() > 0 else None
    chain = atoms.chain.values
    best, pair = -1, None
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            pairs = capped_distance(xyz[chain == a], xyz[chain == b], 4.5, box=box, return_distances=False)
            n = len(np.unique(pairs[:, 1])) if len(pairs) else 0
            if n > best:
                best, pair = n, (a, b)
    return pair


def _sides(u, chains, select, renumber=None):
    """Which atoms make up side a and side b. Returns (atoms, index, (a, b), log)."""
    atoms, index, table = _polymer(u, renumber)
    log = []
    if select:
        if len(select) != 2:
            raise ValueError('Give exactly two selections, e.g. --select "resid 1-110" "resid 111-199"')
        pos = {j: i for i, j in enumerate(index)}
        picked = []
        for label, text in zip('AB', select):
            try:
                found = u.select_atoms(text).indices
            except Exception as e:
                raise ValueError(f'Selection "{text}" not understood: {e}') from e
            rows = [pos[j] for j in found if j in pos]
            if not rows:
                raise ValueError(f'Selection "{text}" holds no protein or nucleic acid atoms.\n{table}')
            picked.append(rows)
        if set(picked[0]) & set(picked[1]):
            raise ValueError('The two selections overlap')
        rows = picked[0] + picked[1]
        atoms, index = atoms.iloc[rows].reset_index(drop=True), index[rows]
        atoms['part'] = atoms.chain          # the chain each atom came from: a side may hold several
        atoms['chain'] = ['A'] * len(picked[0]) + ['B'] * len(picked[1])
        a, b = 'A', 'B'
        log.append(f'Side A = "{select[0]}", side B = "{select[1]}"')
    else:
        present = list(dict.fromkeys(atoms.chain))
        if chains:
            if len(chains) != 2:
                raise ValueError('Give exactly two chains, e.g. --chains A D')
            a, b = chains
            missing = [c for c in (a, b) if c not in present]
            if missing or a == b:
                raise ValueError(f'Chain {" and ".join(missing) or a} not found.\n{table}\n'
                                 'Pick two of these with --chains, or describe the two sides with --select.')
        else:
            a, b = _closest_pair(atoms, u.atoms.positions[index].astype(float), u.trajectory.ts.dimensions)
            log.append(table)
            log.append(f'Following chains {a} and {b} (the pair with most contacts); choose others with --chains')
        if 'lettered' in table and chains:
            log.append(table)
        rows = np.r_[np.where(atoms.chain == a)[0], np.where(atoms.chain == b)[0]]
        atoms, index = atoms.iloc[rows].reset_index(drop=True), index[rows]
    if 'part' in atoms:
        # a side made of several chains (antibody heavy + light, both DNA strands): PandaPPI knows a residue by
        # side and number, so chains whose numbers repeat are shifted by whole thousands, and the shift is said
        for c in (a, b):
            side = atoms.chain == c
            parts = list(dict.fromkeys(atoms.part[side]))
            spans = [set(atoms.resnum[side & (atoms.part == q)]) for q in parts]
            if len(parts) > 1 and sum(map(len, spans)) > len(set().union(*spans)):
                top_, shifts = max(spans[0]), []
                for q, span in zip(parts[1:], spans[1:]):
                    shift = (top_ // 1000 + 1) * 1000 - (min(span) // 1000) * 1000
                    atoms.loc[side & (atoms.part == q), 'resnum'] += shift
                    top_ = max(span) + shift
                    shifts.append(f'{q} +{shift}')
                log.append(f'Side {c} holds {len(parts)} chains with repeating residue numbers: numbers shifted '
                           f'for chain {", ".join(shifts)}')
        atoms = atoms.drop(columns='part')
    if atoms.resnum.max() > 9999 or atoms.resnum.min() < -999:
        raise ValueError('Residue numbers outside -999..9999 cannot be written to the PDB frames PandaPPI reads; '
                         'renumber the topology')
    return atoms.drop(columns='segid'), index, (a, b), log


def _whole(xyz: np.ndarray, b_sel: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Undo periodic wrapping: make each side whole, then put side b at the image closest to side a."""
    inv = np.linalg.inv(box)
    for sel in (~b_sel, b_sel):
        part = xyz[sel]
        # consecutive atoms of a chain are a few Å apart, so any jump of about a box length is a wrap
        part[1:] -= np.cumsum(np.round(np.diff(part, axis=0) @ inv), axis=0) @ box
        xyz[sel] = part
    a, b = xyz[~b_sel][::4], xyz[b_sel][::4]
    base = np.round((b.mean(0) - a.mean(0)) @ inv)
    best, shift = np.inf, base
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                s = base + (dx, dy, dz)
                d = ((a[:, None] - (b - s @ box)[None]) ** 2).sum(-1).min()
                if d < best:
                    best, shift = d, s
    xyz[b_sel] -= shift @ box
    return xyz


def _interface_waters(u, ts, index, atoms, b_sel, xyz, box):
    """Oxygens of waters within hydrogen-bond reach (3.5 Å) of both sides, placed at the image next to the complex."""
    from MDAnalysis.lib.distances import capped_distance
    wat = u._ppitango_water
    if not len(wat):
        return np.zeros((0, 3))
    polar = atoms.element.isin(['N', 'O']).values
    near = []
    for side in (~b_sel, b_sel):
        pairs = capped_distance(ts.positions[index[side & polar]], ts.positions[wat], 3.5, box=ts.dimensions,
                                return_distances=False)
        near.append(set(pairs[:, 1]) if len(pairs) else set())
    both = sorted(near[0] & near[1])
    w = ts.positions[wat[both]].astype(float)
    if len(w) and box is not None:   # put each water at the image closest to the (already whole) complex
        inv = np.linalg.inv(box)
        ref = xyz[polar]
        for k in range(len(w)):
            d = ref - w[k]
            d -= np.round(d @ inv) @ box
            j = (d ** 2).sum(1).argmin()
            w[k] = ref[j] - d[j]
    return w


def _read_trajectory(top: str, traj: Optional[str], chains, stride: int = 1, max_frames: Optional[int] = None,
                     select=None, renumber=None, waters: bool = False):
    """Coordinates of the two sides from a multi-model PDB, or from topology + trajectory (MDAnalysis).

    Returns (atoms, coords, (a, b), times, log, waters); `times` is in ps when the trajectory stores them, else
    None; `waters` is a list of oxygen positions per frame when asked for, else None.
    """
    if traj is None and top.lower().endswith(('.pdb', '.ent')) and not (select or renumber or waters):
        a, b = pick_chains(top, chains)
        atoms, coords, _ = read_pdb(top, (a, b), all_models=True)
        keep = np.arange(len(coords))[::stride]
        if max_frames and len(keep) > max_frames:
            keep = keep[np.linspace(0, len(keep) - 1, max_frames).round().astype(int)]
        order = np.r_[np.where(atoms.chain == a)[0], np.where(atoms.chain == b)[0]]
        return atoms.iloc[order].reset_index(drop=True), coords[keep][:, order], (a, b), None, [], None
    import warnings
    u = _universe(top, traj)
    atoms, index, (a, b), log = _sides(u, chains, select, renumber)
    found = None
    if waters:
        names, resnames = _attr(u.atoms, 'names'), _attr(u.atoms, 'resnames')
        u._ppitango_water = np.where(np.isin(resnames, list(WATER)) & np.isin(names, ['O', 'OW', 'OH2', 'OW1']))[0]
        if not len(u._ppitango_water):
            raise ValueError('--waters: no water molecules in these files (the trajectory was saved without '
                             'solvent?). Give the full-system topology and trajectory.')
        found = []
    keep = np.arange(len(u.trajectory))[::stride]
    if max_frames and len(keep) > max_frames:   # choose the frames first: long trajectories are never read whole
        keep = keep[np.linspace(0, len(keep) - 1, max_frames).round().astype(int)]
    b_sel = (atoms.chain == b).values
    frames, times, unwrapped = [], [], False
    for i in keep:
        ts = u.trajectory[i]
        xyz = ts.positions[index].astype(float)
        box = None
        if ts.dimensions is not None and ts.dimensions[:3].min() > 0:
            before = xyz.copy()
            box = ts.triclinic_dimensions.astype(float)
            xyz = _whole(xyz, b_sel, box)
            unwrapped |= bool(np.abs(xyz - before).max() > 1)
        frames.append(xyz)
        if waters:
            found.append(_interface_waters(u, ts, index, atoms, b_sel, xyz, box))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')   # "Reader has no dt information" for single structures
            times.append(ts.time)
    if unwrapped:
        log.append('Periodic box: chains were made whole and put back next to each other')
    times = np.array(times, dtype=float)
    timed = u.trajectory.format in TIMED and len(times) > 1 and (np.diff(times) > 0).all()
    if waters:
        log.append(f'Interface waters: {np.mean([len(w) for w in found]):.1f} per frame on average sit within '
                   '3.5 Å of both sides')
    return atoms, np.array(frames), (a, b), times if timed else None, log, found


def _fit_atoms(atoms: pd.DataFrame) -> np.ndarray:
    """Atoms frames are superposed on: Cα and nucleotide C4', or everything if there are too few."""
    sel = atoms.name.isin(['CA', "C4'"]).values
    return sel if sel.sum() >= 3 else np.ones(len(atoms), bool)


def md(top: str, traj: Optional[str] = None, chains: Optional[Sequence[str]] = None, stride: int = 1,
       max_frames: Optional[int] = 200, dt: Optional[float] = None, unit: str = 'ns',
       select: Optional[Sequence[str]] = None, renumber: Optional[dict] = None, waters: bool = False) -> Ensemble:
    """A molecular dynamics trajectory.

    `top` is any structure or topology file MDAnalysis reads that holds the same atoms as `traj` (pdb, gro,
    tpr, prmtop/parm7, psf ...); `traj` any trajectory format (xtc, trr, dcd, nc ...). Without `traj`, `top` is
    read as a multi-model PDB. The two sides are `chains` (two chain IDs; without chain IDs in the topology,
    chains are lettered A, B, ... in file order) or `select` (two MDAnalysis selections). `dt` is the time
    between saved frames in `unit`; left out, the times stored in the trajectory are used (xtc, trr, nc), or
    the axis counts frames.
    """
    atoms, coords, (a, b), times, log, wat = _read_trajectory(top, traj, chains, stride, max_frames, select,
                                                              renumber, waters)
    if dt is not None:
        x, xlabel = _frame_numbers(top, traj, stride, max_frames, len(coords)) * dt, f'Time ({unit})'
    elif times is not None:
        unit = 'ns' if times[-1] - times[0] >= 1000 else 'ps'
        x, xlabel = (times - times[0]) * (1e-3 if unit == 'ns' else 1), f'Time ({unit})'
    else:
        x, xlabel, unit = _frame_numbers(top, traj, stride, max_frames, len(coords)).astype(float), 'Frame', 'frames'
        log.append('No frame times known: the axis counts frames. Give --dt (time between saved frames) for a '
                   'time axis')
    coords = _superpose(coords, _fit_atoms(atoms), wat)
    return Ensemble(atoms, coords, (a, b), x=x, xlabel=xlabel, kind='md', waters=wat, log=log,
                    title=f'{_name(traj or top)} · molecular dynamics',
                    notes=[f'{len(coords)} frames over {x[-1] - x[0]:g} {unit}'])


def _frame_numbers(top, traj, stride, max_frames, n_kept) -> np.ndarray:
    """Index in the trajectory of every kept frame (the same choice _read_trajectory makes)."""
    if traj is None and top.lower().endswith(('.pdb', '.ent')):
        with open(top) as f:
            n = sum(line.startswith('ENDMDL') for line in f) or 1
    else:
        n = len(_universe(top, traj).trajectory)
    keep = np.arange(n)[::stride]
    if max_frames and len(keep) > max_frames:
        keep = keep[np.linspace(0, len(keep) - 1, max_frames).round().astype(int)]
    assert len(keep) == n_kept
    return keep


def fep(top: str, windows: Sequence[str], chains: Optional[Sequence[str]] = None, dg: Optional[str] = None,
        frames_per_window: int = 8, label: str = '', select: Optional[Sequence[str]] = None,
        renumber: Optional[dict] = None, waters: bool = False) -> Ensemble:
    """Free energy perturbation: one trajectory (or multi-model PDB) per lambda window, in order.

    `dg` is a CSV with one row per window and a `dG` column (cumulative free energy, kcal/mol);
    an optional `lambda` column sets the axis, otherwise windows are spaced evenly from 0 to 1.
    """
    stacks, lam, log, wat = [], [], [], [] if waters else None
    table = pd.read_csv(dg) if dg else None
    if table is not None and len(table) != len(windows):
        raise ValueError(f'{dg} has {len(table)} rows but there are {len(windows)} windows')
    if table is not None and 'dG' not in table:
        raise ValueError(f'{dg} needs a column called dG (cumulative free energy per window); it has '
                         f'{list(table.columns)}')
    lambdas = table['lambda'].values if table is not None and 'lambda' in table else \
        np.linspace(0, 1, len(windows))
    for k, w in enumerate(windows):
        is_pdb = w.lower().endswith(('.pdb', '.ent'))
        atoms, coords, pair, _, log_k, wat_k = _read_trajectory(w if is_pdb else top, None if is_pdb else w,
                                                                chains, 1, frames_per_window, select, renumber,
                                                                waters)
        if waters:
            wat += wat_k
        chains = chains or (None if select else pair)   # later windows follow the pair found in the first
        log = log or log_k
        stacks.append(coords)
        # frames of a window fill the stretch up to the next lambda, so the axis keeps moving;
        # the last window stays at the final lambda
        nxt = lambdas[k + 1] if k + 1 < len(lambdas) else lambdas[k]
        lam.append(lambdas[k] + (nxt - lambdas[k]) * np.arange(len(coords)) / len(coords))
    coords = _superpose(np.concatenate(stacks), _fit_atoms(atoms), wat)
    x = np.concatenate(lam)
    curve = np.interp(x, lambdas, table['dG'].values) if table is not None else None
    return Ensemble(atoms, coords, pair, x=x, xlabel='λ (alchemical progress)', kind='fep', waters=wat, log=log,
                    title=f'{_name(top)} · free energy perturbation' + (f' · {label}' if label else ''),
                    curve=curve, curve_label='ΔΔG (kcal/mol)' if table is not None else '',
                    notes=[f'{len(windows)} λ windows, {len(coords)} frames',
                           'contacts are geometric: decoupled atoms still count'])
