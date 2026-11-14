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
    log: List[str] = field(default_factory=list)     # what was decided while reading (chains found, box handling)

    @property
    def n_frames(self) -> int:
        return len(self.coords)

    def pdb(self, frame: int) -> str:
        """One frame as PDB text."""
        return _pdb_text(self.atoms, self.coords[frame])


def _pdb_text(atoms: pd.DataFrame, xyz: np.ndarray) -> str:
    lines = []
    for i, (a, (x, y, z)) in enumerate(zip(atoms.itertuples(), xyz), 1):
        name = a.name if len(a.name) == 4 else f' {a.name:<3}'
        lines.append(f'ATOM  {i % 100000:5d} {name} {a.resname:>3} {a.chain}{a.resnum:4d}    '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a.element:>2}')
    return '\n'.join(lines) + '\nEND\n'


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
                rows.append((chain, RESNAME_ALIASES.get(resname, resname), int(line[22:26]), name, element))
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


def _superpose(coords: np.ndarray, sel: np.ndarray) -> np.ndarray:
    """Fit every frame onto the first, using the selected atoms (Kabsch)."""
    ref = coords[0, sel]
    ref_c = ref.mean(0)
    out = np.empty_like(coords)
    for i, xyz in enumerate(coords):
        c = xyz[sel].mean(0)
        u, _, vt = np.linalg.svd((xyz[sel] - c).T @ (ref - ref_c))
        d = np.sign(np.linalg.det(u @ vt))
        rot = u @ np.diag([1, 1, d]) @ vt
        out[i] = (xyz - c) @ rot + ref_c
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
        frames: int = 48, cutoff: float = 15.0) -> Ensemble:
    """One oscillation along an elastic-network normal mode of the two-chain complex.

    Each residue moves rigidly with its Cα, and the mode is scaled so the Cα RMSD from the
    starting structure peaks at `rmsd` Å.
    """
    a, b = pick_chains(pdb, chains)
    atoms, coords, _ = read_pdb(pdb, (a, b))
    xyz = coords[0]
    res_key = atoms.chain + atoms.resnum.astype(str)
    is_ca = (atoms.name == 'CA').values
    ca_of = dict(zip(res_key[is_ca], range(is_ca.sum())))
    missing = ~res_key.isin(ca_of)
    if missing.any():  # residues without a Cα (nucleotides, truncated residues) cannot be moved
        atoms, xyz, res_key, is_ca = atoms[~missing.values].reset_index(drop=True), xyz[~missing.values], \
            res_key[~missing.values], is_ca[~missing.values]
    gone = [c for c in (a, b) if not (atoms.chain == c).any()]
    if gone:
        raise ValueError(f'Normal modes move residues with their Cα atom and chain {" and ".join(gone)} has none '
                         '(a nucleic acid?). Choose two protein chains with --chains.')
    vals, modes = anm_modes(xyz[is_ca], cutoff, n_modes=max(mode, 1))
    disp = modes[mode - 1]
    disp = disp / np.sqrt((disp ** 2).sum(-1).mean()) * rmsd      # Cα RMSD = rmsd at full swing
    per_atom = disp[res_key.map(ca_of).values]
    phase = np.arange(frames) / frames
    stack = xyz[None] + np.sin(2 * np.pi * phase)[:, None, None] * per_atom[None]
    return Ensemble(atoms, stack, (a, b), x=phase, xlabel=f'Phase of normal mode {mode} (cycles)', kind='nma',
                    title=f'{_name(pdb)} · normal mode {mode}', periodic=True,
                    notes=[f'Elastic network model (Cα, {cutoff:g} Å cutoff), mode {mode}, ±{rmsd:g} Å Cα RMSD'])


# ---------------------------------------------------------------- trajectories

# force-field atom names -> the PDB names PandaPPI's detectors know
ATOM_ALIASES = {('ILE', 'CD'): 'CD1', 'OT1': 'O', 'OC1': 'O', 'O1': 'O', 'OT2': 'OXT', 'OC2': 'OXT', 'O2': 'OXT'}
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


def _polymer(u):
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

    resname = _attr(ag, 'resnames')[index]
    name = names[index]
    resname = np.array([RESNAME_ALIASES.get(r, r) for r in resname])
    name = np.array([ATOM_ALIASES.get((r, n), n) for r, n in zip(resname, name)])
    # a C-terminal oxygen pair is O/OXT only on protein residues; nucleic acids keep O2 and friends
    is_prot = protein[resix[index]]
    name = np.where(is_prot, [ATOM_ALIASES.get(n, n) for n in name], name)
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


def _sides(u, chains, select):
    """Which atoms make up side a and side b. Returns (atoms, index, (a, b), log)."""
    atoms, index, table = _polymer(u)
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


def _read_trajectory(top: str, traj: Optional[str], chains, stride: int = 1, max_frames: Optional[int] = None,
                     select=None):
    """Coordinates of the two sides from a multi-model PDB, or from topology + trajectory (MDAnalysis).

    Returns (atoms, coords, (a, b), times, log); `times` is in ps when the trajectory stores them, else None.
    """
    if traj is None and top.lower().endswith(('.pdb', '.ent')) and not select:
        a, b = pick_chains(top, chains)
        atoms, coords, _ = read_pdb(top, (a, b), all_models=True)
        keep = np.arange(len(coords))[::stride]
        if max_frames and len(keep) > max_frames:
            keep = keep[np.linspace(0, len(keep) - 1, max_frames).round().astype(int)]
        order = np.r_[np.where(atoms.chain == a)[0], np.where(atoms.chain == b)[0]]
        return atoms.iloc[order].reset_index(drop=True), coords[keep][:, order], (a, b), None, []
    import warnings
    u = _universe(top, traj)
    atoms, index, (a, b), log = _sides(u, chains, select)
    keep = np.arange(len(u.trajectory))[::stride]
    if max_frames and len(keep) > max_frames:   # choose the frames first: long trajectories are never read whole
        keep = keep[np.linspace(0, len(keep) - 1, max_frames).round().astype(int)]
    b_sel = (atoms.chain == b).values
    frames, times, unwrapped = [], [], False
    for i in keep:
        ts = u.trajectory[i]
        xyz = ts.positions[index].astype(float)
        if ts.dimensions is not None and ts.dimensions[:3].min() > 0:
            before = xyz.copy()
            xyz = _whole(xyz, b_sel, ts.triclinic_dimensions.astype(float))
            unwrapped |= bool(np.abs(xyz - before).max() > 1)
        frames.append(xyz)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')   # "Reader has no dt information" for single structures
            times.append(ts.time)
    if unwrapped:
        log.append('Periodic box: chains were made whole and put back next to each other')
    times = np.array(times, dtype=float)
    timed = u.trajectory.format in TIMED and len(times) > 1 and (np.diff(times) > 0).all()
    return atoms, np.array(frames), (a, b), times if timed else None, log


def _fit_atoms(atoms: pd.DataFrame) -> np.ndarray:
    """Atoms frames are superposed on: Cα and nucleic acid P, or everything if there are too few."""
    sel = atoms.name.isin(['CA', 'P']).values
    return sel if sel.sum() >= 3 else np.ones(len(atoms), bool)


def md(top: str, traj: Optional[str] = None, chains: Optional[Sequence[str]] = None, stride: int = 1,
       max_frames: Optional[int] = 200, dt: Optional[float] = None, unit: str = 'ns',
       select: Optional[Sequence[str]] = None) -> Ensemble:
    """A molecular dynamics trajectory.

    `top` is any structure or topology file MDAnalysis reads that holds the same atoms as `traj` (pdb, gro,
    tpr, prmtop/parm7, psf ...); `traj` any trajectory format (xtc, trr, dcd, nc ...). Without `traj`, `top` is
    read as a multi-model PDB. The two sides are `chains` (two chain IDs; without chain IDs in the topology,
    chains are lettered A, B, ... in file order) or `select` (two MDAnalysis selections). `dt` is the time
    between saved frames in `unit`; left out, the times stored in the trajectory are used (xtc, trr, nc), or
    the axis counts frames.
    """
    atoms, coords, (a, b), times, log = _read_trajectory(top, traj, chains, stride, max_frames, select)
    if dt is not None:
        x, xlabel = _frame_numbers(top, traj, stride, max_frames, len(coords)) * dt, f'Time ({unit})'
    elif times is not None:
        unit = 'ns' if times[-1] - times[0] >= 1000 else 'ps'
        x, xlabel = (times - times[0]) * (1e-3 if unit == 'ns' else 1), f'Time ({unit})'
    else:
        x, xlabel, unit = _frame_numbers(top, traj, stride, max_frames, len(coords)).astype(float), 'Frame', 'frames'
        log.append('No frame times known: the axis counts frames. Give --dt (time between saved frames) for a '
                   'time axis')
    coords = _superpose(coords, _fit_atoms(atoms))
    ens = Ensemble(atoms, coords, (a, b), x=x, xlabel=xlabel, kind='md',
                   title=f'{_name(traj or top)} · molecular dynamics',
                   notes=[f'{len(coords)} frames over {x[-1] - x[0]:g} {unit}'])
    ens.log = log
    return ens


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
        frames_per_window: int = 8, label: str = '', select: Optional[Sequence[str]] = None) -> Ensemble:
    """Free energy perturbation: one trajectory (or multi-model PDB) per lambda window, in order.

    `dg` is a CSV with one row per window and a `dG` column (cumulative free energy, kcal/mol);
    an optional `lambda` column sets the axis, otherwise windows are spaced evenly from 0 to 1.
    """
    stacks, lam, log = [], [], []
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
        atoms, coords, pair, _, log_k = _read_trajectory(w if is_pdb else top, None if is_pdb else w, chains, 1,
                                                         frames_per_window, select)
        chains = chains or (None if select else pair)   # later windows follow the pair found in the first
        log = log or log_k
        stacks.append(coords)
        # frames of a window fill the stretch up to the next lambda, so the axis keeps moving;
        # the last window stays at the final lambda
        nxt = lambdas[k + 1] if k + 1 < len(lambdas) else lambdas[k]
        lam.append(lambdas[k] + (nxt - lambdas[k]) * np.arange(len(coords)) / len(coords))
    coords = _superpose(np.concatenate(stacks), _fit_atoms(atoms))
    x = np.concatenate(lam)
    curve = np.interp(x, lambdas, table['dG'].values) if table is not None else None
    ens = Ensemble(atoms, coords, pair, x=x, xlabel='λ (alchemical progress)', kind='fep',
                   title=f'{_name(top)} · free energy perturbation' + (f' · {label}' if label else ''),
                   curve=curve, curve_label='ΔΔG (kcal/mol)' if table is not None else '',
                   notes=[f'{len(windows)} λ windows, {len(coords)} frames',
                          'contacts are geometric: decoupled atoms still count'])
    ens.log = log
    return ens
