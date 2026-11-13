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

def _read_trajectory(top: str, traj: Optional[str], chains, stride: int, max_frames: Optional[int]):
    """Coordinates of the two chains from a multi-model PDB, or from topology + trajectory (MDAnalysis)."""
    if traj is None:
        atoms, coords, _ = read_pdb(top, chains, all_models=True)
        coords = coords[::stride]
    else:
        try:
            import MDAnalysis as mda
        except ImportError as e:
            raise ImportError('Reading trajectory files needs MDAnalysis: pip install "ppi-tango[md]"') from e
        import warnings
        atoms, _, index = read_pdb(top, chains)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')  # MDAnalysis is chatty about PDB records and future changes
            u = mda.Universe(top, traj)
        b_sel = (atoms.chain == chains[1]).values
        frames = []
        for ts in u.trajectory[::stride]:
            xyz = ts.positions[index].astype(float)
            if ts.dimensions is not None and ts.dimensions[:3].min() > 0:
                # bring chain b back next to chain a if the periodic box split the complex
                box = ts.triclinic_dimensions.astype(float)
                gap = xyz[b_sel].mean(0) - xyz[~b_sel].mean(0)
                xyz[b_sel] -= np.round(gap @ np.linalg.inv(box)) @ box
            frames.append(xyz)
        coords = np.array(frames)
    if max_frames and len(coords) > max_frames:
        coords = coords[np.linspace(0, len(coords) - 1, max_frames).round().astype(int)]
    return atoms, coords


def md(top: str, traj: Optional[str] = None, chains: Optional[Sequence[str]] = None, stride: int = 1,
       max_frames: Optional[int] = 200, dt: float = 1.0, unit: str = 'ns') -> Ensemble:
    """A molecular dynamics trajectory. `top` is a PDB with chain IDs; `traj` any format MDAnalysis
    reads (xtc, dcd, trr, nc ...). Without `traj`, `top` is read as a multi-model PDB.
    `dt` is the time between saved frames, in `unit`."""
    a, b = pick_chains(top, chains)
    atoms, coords = _read_trajectory(top, traj, (a, b), stride, None)
    n = len(coords)
    x = np.arange(n) * dt * stride
    if max_frames and n > max_frames:
        keep = np.linspace(0, n - 1, max_frames).round().astype(int)
        coords, x = coords[keep], x[keep]
    coords = _superpose(coords, (atoms.name == 'CA').values)
    return Ensemble(atoms, coords, (a, b), x=x, xlabel=f'Time ({unit})', kind='md',
                    title=f'{_name(traj or top)} · molecular dynamics',
                    notes=[f'{len(coords)} frames over {x[-1]:g} {unit}'])


def fep(top: str, windows: Sequence[str], chains: Optional[Sequence[str]] = None, dg: Optional[str] = None,
        frames_per_window: int = 8, label: str = '') -> Ensemble:
    """Free energy perturbation: one trajectory (or multi-model PDB) per lambda window, in order.

    `dg` is a CSV with one row per window and a `dG` column (cumulative free energy, kcal/mol);
    an optional `lambda` column sets the axis, otherwise windows are spaced evenly from 0 to 1.
    """
    a, b = pick_chains(top, chains)
    stacks, lam = [], []
    table = pd.read_csv(dg) if dg else None
    if table is not None and len(table) != len(windows):
        raise ValueError(f'{dg} has {len(table)} rows but there are {len(windows)} windows')
    lambdas = table['lambda'].values if table is not None and 'lambda' in table else \
        np.linspace(0, 1, len(windows))
    for k, w in enumerate(windows):
        is_pdb = w.lower().endswith(('.pdb', '.ent'))
        atoms, coords = _read_trajectory(w if is_pdb else top, None if is_pdb else w, (a, b), 1, frames_per_window)
        stacks.append(coords)
        # frames of a window fill the stretch up to the next lambda, so the axis keeps moving;
        # the last window stays at the final lambda
        nxt = lambdas[k + 1] if k + 1 < len(lambdas) else lambdas[k]
        lam.append(lambdas[k] + (nxt - lambdas[k]) * np.arange(len(coords)) / len(coords))
    coords = _superpose(np.concatenate(stacks), (atoms.name == 'CA').values)
    x = np.concatenate(lam)
    curve = np.interp(x, lambdas, table['dG'].values) if table is not None else None
    return Ensemble(atoms, coords, (a, b), x=x, xlabel='λ (alchemical progress)', kind='fep',
                    title=f'{_name(top)} · free energy perturbation' + (f' · {label}' if label else ''),
                    curve=curve, curve_label='ΔΔG (kcal/mol)' if table is not None else '',
                    notes=[f'{len(windows)} λ windows, {len(coords)} frames',
                           'contacts are geometric: decoupled atoms still count'])
