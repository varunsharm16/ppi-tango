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
    sides: tuple = ()              # what the two sides are called: ('Chain A', 'Chain D'), ('Chain C', 'Chains A+B')

    def __post_init__(self):
        if not self.sides:
            self.sides = tuple(f'Chain {c}' for c in self.chains)
        if 'resid' not in self.atoms:      # the name a residue is shown by; `resnum` is what the PDB frames carry
            self.atoms = self.atoms.assign(resid=self.atoms.resnum.astype(str))
        if 'kind' not in self.atoms:
            self.atoms = self.atoms.assign(kind='protein')

    @property
    def labels(self) -> pd.Series:
        """Every atom's residue as it is shown: 'TYR 100A', or 'TYR H:100A' when a side holds several chains."""
        return self.atoms.resname + ' ' + self.atoms.resid

    @property
    def n_frames(self) -> int:
        return len(self.coords)

    @property
    def pair(self) -> str:
        """The two sides in a sentence: 'chains A and D', 'chain C and chains A+B'."""
        a, b = self.sides
        if a.startswith('Chain ') and b.startswith('Chain '):
            return f'chains {a[6:]} and {b[6:]}'
        lower = lambda s: s[0].lower() + s[1:] if s.startswith('Chain') else s
        return f'{lower(a)} and {lower(b)}'

    def pdb(self, frame: int) -> str:
        """One frame as PDB text."""
        return _pdb_text(self.atoms, self.coords[frame], self.waters[frame] if self.waters is not None else None)


def _pdb_text(atoms: pd.DataFrame, xyz: np.ndarray, waters: Optional[np.ndarray] = None) -> str:
    lines = []
    for i, (a, (x, y, z)) in enumerate(zip(atoms.itertuples(), xyz), 1):
        name = a.name if len(a.name) == 4 else f' {a.name:<3}'
        lines.append(f'ATOM  {i % 100000:5d} {name[:4]} {a.resname[:3]:>3} {a.chain}{a.resnum:4d}    '
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


def _number(atoms: pd.DataFrame) -> pd.DataFrame:
    """Settle how residues are numbered. `resnum` is the number the PDB frames carry: it has to be unique within
    a side and fit four digits, so a side with insertion codes (100, 100A), with chains whose numbers repeat, or
    with numbers above 9999 is counted 1, 2, 3 ... instead. `resid` is what people read and never changes: the
    original number with its insertion code, and the chain in front when a side holds several (H:100A)."""
    atoms = atoms.copy()
    icode = atoms.icode if 'icode' in atoms else pd.Series('', index=atoms.index)
    part = atoms.part if 'part' in atoms else atoms.chain
    atoms['resid'] = atoms.resnum.astype(str) + icode
    for c in dict.fromkeys(atoms.chain):
        side = (atoms.chain == c).values
        key = pd.DataFrame({'p': part[side], 'n': atoms.resnum[side], 'i': icode[side]})
        new = (key != key.shift()).any(axis=1).values            # first atom of every residue
        if part[side].nunique() > 1:
            atoms.loc[side, 'resid'] = part[side] + ':' + atoms.resid[side]
        nums = atoms.resnum[side].values
        if len(set(nums[new])) < new.sum() or nums.max() > 9999 or nums.min() < -999:
            if new.sum() > 9999:
                raise ValueError(f'Side {c} holds {new.sum()} residues; the PDB frames PandaPPI reads take 9999. '
                                 'Choose a smaller part with --select')
            atoms.loc[side, 'resnum'] = new.cumsum()
    return atoms.drop(columns=[c for c in ('icode', 'part') if c in atoms])


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
                rows.append((chain, resname, int(line[22:26]), name, element, line[26].strip()))
                index.append(n)
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if xyz:
        models.append(xyz)
    atoms = pd.DataFrame(rows, columns=['chain', 'resname', 'resnum', 'name', 'element', 'icode'])
    if atoms.empty:
        raise ValueError(f'No polymer atoms found in {path}' + (f' for chains {list(chains)}' if chains else ''))
    atoms = _number(atoms)
    if len({len(m) for m in models}) != 1:
        raise ValueError(f'{path}: models do not all contain the same atoms')
    return atoms, np.array(models, dtype=float), np.array(index)


def pick_chains(path: str, chains: Optional[Sequence[str]]) -> tuple:
    """The two chains to follow: the ones asked for, or the pair with the most atoms in contact."""
    if chains:
        if len(chains) != 2 or chains[0] == chains[1]:
            raise ValueError('Give two different chains, e.g. --chains A D')
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


def _plain_pdb(path: str) -> bool:
    """A PDB file the built-in reader is enough for: protein chains only (nucleic acids go through MDAnalysis,
    which sorts chains into protein and nucleic acid sides)."""
    if not path.lower().endswith(('.pdb', '.ent')):
        return False
    with open(path, errors='replace') as f:
        return not any(line.startswith('ATOM') and line[12:16].strip() in ("C4'", 'C4*') for line in f)


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
    log, sides = [], ()
    if mode < 1:
        raise ValueError(f'--mode counts from 1 (the slowest mode); got {mode}')
    if select or renumber or not (_plain_pdb(pdb) or (chains and pdb.lower().endswith(('.pdb', '.ent')))):
        u = _universe(pdb, None)
        atoms, index, (a, b), log, sides = _sides(u, chains, select, renumber)
        xyz = u.atoms.positions[index].astype(float)
        box = u.trajectory.ts.dimensions
        if box is not None and box[:3].min() > 0:
            xyz = _whole(xyz, atoms.unit.values, u.trajectory.ts.triclinic_dimensions.astype(float))
    else:
        a, b = pick_chains(pdb, chains)
        atoms, coords, _ = read_pdb(pdb, (a, b))
        xyz = coords[0]
    res_key = atoms.chain + atoms.resnum.astype(str)
    is_node = atoms.name.isin(['CA', "C4'"]).values
    if 'kind' in atoms:   # a ligand, ion or sugar has no backbone: its first atom is its node
        is_node = is_node | ((atoms.kind == 'other') & ~res_key.duplicated()).values
    node_of = dict(zip(res_key[is_node], range(is_node.sum())))
    missing = ~res_key.isin(node_of)
    if missing.any():  # residues without a node atom (truncated residues) cannot be moved
        atoms, xyz, res_key, is_node = atoms[~missing.values].reset_index(drop=True), xyz[~missing.values], \
            res_key[~missing.values], is_node[~missing.values]
    gone = [c for c in (a, b) if not (atoms.chain == c).any()]
    if gone:
        raise ValueError(f"Chain {' and '.join(gone)} has no Cα or C4' atoms to build the elastic network on")
    if mode > 3 * is_node.sum() - 6:
        raise ValueError(f'--mode {mode}: the elastic network of {is_node.sum()} residues has '
                         f'{3 * is_node.sum() - 6} modes')
    vals, modes = anm_modes(xyz[is_node], cutoff, n_modes=mode)
    disp = modes[mode - 1]
    disp = disp / np.sqrt((disp ** 2).sum(-1).mean()) * rmsd      # node RMSD = rmsd at full swing
    per_atom = disp[res_key.map(node_of).values]
    phase = np.arange(frames) / frames
    stack = xyz[None] + np.sin(2 * np.pi * phase)[:, None, None] * per_atom[None]
    nodes = "Cα" + (" and C4'" if (atoms.name == "C4'").any() else '')
    return Ensemble(atoms, stack, (a, b), x=phase, xlabel=f'Phase of normal mode {mode} (cycles)', kind='nma',
                    title=f'{_name(pdb)} · normal mode {mode}', periodic=True, log=log, sides=sides,
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


def _gromacs_include_dir():
    """The force field folder of a GROMACS installation (share/gromacs/top), which a .top includes files from."""
    import os
    import shutil
    found = [os.environ.get('GMXLIB'), os.path.join(os.environ.get('GMXDATA', ''), 'top')]
    for exe in ('gmx', 'gmx_mpi', 'gmx_d'):
        path = shutil.which(exe)
        if path:
            for start in (path, os.path.realpath(path)):
                found.append(os.path.join(os.path.dirname(os.path.dirname(start)), 'share', 'gromacs', 'top'))
    return next((d for d in found if d and os.path.isdir(d)), None)


def _is_gromacs_top(path: str) -> bool:
    """.top is the extension of both a GROMACS topology and an Amber prmtop; the Amber one starts with %VERSION."""
    if not path.lower().endswith(('.top', '.itp')):
        return False
    with open(path, errors='replace') as f:
        return not f.read(200).lstrip().startswith('%')


IONS = {'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'LI', 'RB', 'CS', 'NA+', 'CL-', 'K+', 'SOD', 'CLA', 'POT', 'CAL', 'MN',
        'FE', 'CU', 'NI', 'CO', 'CD', 'BR', 'I', 'F'}
WATER_SITES = {3: ('OW', 'HW1', 'HW2'), 4: ('OW', 'HW1', 'HW2', 'MW'), 5: ('OW', 'HW1', 'HW2', 'LP1', 'LP2')}


def _gromacs_top(path: str, n_atoms: Optional[int] = None):
    """Atoms of a GROMACS .top, in the order of the system: (names, resnames, resnums, molecule of each atom, note).

    Reads the [ moleculetype ], [ atoms ] and [ molecules ] sections and follows #include, #ifdef and #define.
    Residue numbers are those of each molecule (GROMACS itself counts on through the system). A molecule defined
    in a force field file that is not at hand (water, ions: no GROMACS on this machine) is filled in when it can
    be nothing else: ions hold one atom, and the water model is the one that makes the trajectory's `n_atoms`.
    """
    import os
    include = _gromacs_include_dir()
    types, order, missing = {}, [], []
    state = {'section': '', 'mol': None}

    def read(file, defined, depth=0):
        folder = os.path.dirname(os.path.abspath(file))
        active = [True]
        with open(file, errors='replace') as f:
            lines = f.read().replace('\\\n', ' ').splitlines()
        for line in lines:
            line = line.split(';')[0].strip()
            if not line:
                continue
            if line.startswith('#'):
                word, _, rest = line[1:].strip().partition(' ')
                rest = rest.strip()
                if word in ('ifdef', 'ifndef'):
                    active.append(active[-1] and ((rest.split()[0] in defined) == (word == 'ifdef')))
                elif word == 'else':
                    active[-1] = not active[-1] and all(active[:-1])
                elif word == 'endif':
                    active.pop()
                elif not active[-1]:
                    pass
                elif word == 'define' and rest:
                    defined.add(rest.split()[0])
                elif word == 'include' and depth < 20:
                    name = rest.strip('"<>')
                    for base in (folder, os.path.dirname(os.path.abspath(path)), include):
                        if base and os.path.isfile(os.path.join(base, name)):
                            read(os.path.join(base, name), defined, depth + 1)
                            break
                    else:
                        missing.append(name)
                continue
            if not active[-1]:
                continue
            if line.startswith('['):
                state['section'] = line.strip('[] ').lower()
                continue
            cols = line.split()
            if state['section'] == 'moleculetype':
                state['mol'] = cols[0]
                types[cols[0]] = []
            elif state['section'] == 'atoms' and state['mol'] and len(cols) >= 5:
                types[state['mol']].append((cols[4], cols[3], int(cols[2])))      # atom, residue, residue number
            elif state['section'] == 'molecules' and len(cols) >= 2:
                order.append((cols[0], int(cols[1])))

    read(path, set())
    if not order:
        raise ValueError('no [ molecules ] section')
    unknown = [m for m, _ in order if m not in types]
    note = ''
    if unknown:
        waters = sum(n for m, n in order if m not in types and m.upper() in WATER)
        for m in set(unknown):
            if m.upper() in IONS:
                types[m] = [(m, m, 1)]
            elif m.upper() not in WATER:
                raise ValueError(f'molecule {m} is defined in a file that was not found ({", ".join(missing)})')
        if waters:
            known = sum(len(types[m]) * n for m, n in order if m in types)
            sites = (n_atoms - known) / waters if n_atoms else 3
            if n_atoms and n_atoms <= known:      # a trajectory saved without its water: the model does not matter
                sites = 3
            if sites not in WATER_SITES:
                raise ValueError(f'cannot tell the water model: the force field files ({", ".join(missing)}) '
                                 'were not found and the atom count does not settle it')
            for m in set(unknown):
                if m.upper() in WATER:
                    types[m] = [(name, 'SOL', 1) for name in WATER_SITES[int(sites)]]
        note = (f'{os.path.basename(path)}: force field files not found ({", ".join(dict.fromkeys(missing))}); '
                'water and ions were filled in from their names')
    names, resnames, resnums, mol = [], [], [], []
    k = 0
    for m, n in order:
        atoms = types[m]
        for _ in range(n):
            names += [a[0] for a in atoms]
            resnames += [a[1] for a in atoms]
            resnums += [a[2] for a in atoms]
            mol += [k] * len(atoms)
            k += 1
    return names, resnames, np.array(resnums), np.array(mol), [m for m, n in order for _ in range(n)], note


def _empty(mda, names, resnames, resnums, groups, segids=None, chain_ids=None, icodes=None):
    """A Universe from plain per-atom arrays; a new residue starts where `groups` or the residue number changes."""
    resnums, groups = np.asarray(resnums), np.asarray(groups)
    start = np.r_[True, (resnums[1:] != resnums[:-1]) | (groups[1:] != groups[:-1])]
    resindex = np.cumsum(start) - 1
    first = np.where(start)[0]
    seg_of_res = groups[first]
    seg_start = np.r_[True, seg_of_res[1:] != seg_of_res[:-1]]
    u = mda.Universe.empty(len(names), n_residues=len(first), n_segments=int(seg_start.sum()),
                           atom_resindex=resindex, residue_segindex=np.cumsum(seg_start) - 1, trajectory=False)
    u.add_TopologyAttr('names', list(names))
    u.add_TopologyAttr('resnames', [resnames[i] for i in first])
    u.add_TopologyAttr('resids', resnums[first])
    if segids is not None:
        u.add_TopologyAttr('segids', [segids[i] for i in first[seg_start]])
    if chain_ids is not None:
        u.add_TopologyAttr('chainIDs', list(chain_ids))
    if icodes is not None:
        u.add_TopologyAttr('icodes', [icodes[i] for i in first])
    return u


def _n_atoms(mda, traj: str) -> int:
    from MDAnalysis.coordinates.core import reader
    r = reader(traj)
    try:
        return r.n_atoms
    finally:
        r.close()


def _topology(mda, top: str, n_traj: Optional[int]):
    """The Universe of a topology file, without the trajectory, and what had to be said about reading it."""
    import os
    notes = []
    if _is_gromacs_top(top):
        try:
            names, resnames, resnums, mol, molname, note = _gromacs_top(top, n_traj)
            notes += [note] if note else []
            segids = np.array(molname)[mol]
            return _empty(mda, names, resnames, resnums, mol, segids=segids), notes
        except (ValueError, IndexError, OSError) as own:     # an unusual .top: let MDAnalysis try
            include = _gromacs_include_dir()
            try:
                return mda.Universe(top, topology_format='ITP', **({'include_dir': include} if include else {})), notes
            except Exception:
                raise ValueError(f'{own}. A GROMACS .top includes force field files from the GROMACS '
                                 'installation: load GROMACS first (gmx on the PATH, or GMXDATA set), keep the '
                                 '.itp files it includes next to it, or give the .gro or .tpr instead') from own
    try:
        return mda.Universe(top), notes
    except Exception as e:
        if 'tpx version' not in str(e):
            raise
        # MDAnalysis lags behind new GROMACS releases. The .gro mdrun writes next to the .tpr holds the same
        # atoms and never changes format, and so does the .top the .tpr was made from
        name, stem, folder = os.path.basename(top), os.path.splitext(top)[0], os.path.dirname(top)
        for other in (stem + '.gro', stem + '.top', os.path.join(folder, 'topol.top')):
            if os.path.exists(other):
                try:
                    u, more = _topology(mda, other, n_traj)
                except Exception:
                    continue
                return u, [f'{name} is from a GROMACS newer than the installed MDAnalysis reads: used '
                           f'{os.path.basename(other)} next to it instead'] + more
        raise ValueError(f'{e}. This .tpr is from a GROMACS newer than the installed MDAnalysis reads. Give the '
                         f'.gro or the .top of the same system instead (gmx editconf -f {name} -o system.gro '
                         'writes a .gro), or try pip install -U MDAnalysis') from e


def _original_numbers(mda, u, top: str, notes: list):
    """Residue numbers as the structure had them, where the topology format dropped them.

    A .tpr counts residues 1, 2, 3 ... through the whole system: the .gro of the same run next to it still has
    the numbers of each chain. An Amber prmtop does the same unless it carries the original numbers and chain
    IDs (written by parmed's addPDB, or tleap from a PDB), which MDAnalysis does not read."""
    import os
    low = top.lower()
    if low.endswith('.tpr'):
        gro = os.path.splitext(top)[0] + '.gro'
        if os.path.exists(gro):
            try:
                g = mda.Universe(gro)
                if len(g.atoms) == len(u.atoms) and len(g.residues) == len(u.residues):
                    u.residues.resids = g.residues.resids
                    notes.append(f'Residue numbers taken from {os.path.basename(gro)} (a .tpr counts straight '
                                 'through all chains)')
            except Exception:
                pass
    elif low.endswith(('.prmtop', '.parm7', '.top')) and not _is_gromacs_top(top):
        flags, flag = {}, None
        with open(top, errors='replace') as f:
            for line in f:
                if line.startswith('%FLAG'):
                    flag = line.split()[1]
                    flag = flag if flag in ('RESIDUE_NUMBER', 'RESIDUE_CHAINID', 'RESIDUE_ICODE') else None
                    if flag:
                        flags[flag] = []
                elif flag and line.startswith('%FORMAT'):
                    import re
                    flags[flag].append(int(re.search(r'\d+[aAiI](\d+)', line).group(1)))
                elif flag:
                    w = flags[flag][0]
                    row = line.rstrip('\n')
                    flags[flag] += [row[i:i + w] for i in range(0, len(row), w)]
        n = len(u.residues)
        got = {k: [x.strip() for x in v[1:]] for k, v in flags.items() if len(v) - 1 == n}
        if 'RESIDUE_NUMBER' in got:
            u.residues.resids = np.array(got['RESIDUE_NUMBER'], dtype=int)
            notes.append('Residue numbers taken from the RESIDUE_NUMBER section of the prmtop')
        if 'RESIDUE_CHAINID' in got and len(set(got['RESIDUE_CHAINID']) - {''}) >= 2:
            u.add_TopologyAttr('chainIDs', np.array(got['RESIDUE_CHAINID'])[u.atoms.resindices])
        if 'RESIDUE_ICODE' in got and any(got['RESIDUE_ICODE']):
            u.add_TopologyAttr('icodes', got['RESIDUE_ICODE'])


def _groups(u) -> list:
    """Sets of atoms a trajectory is commonly saved for, smallest first: (description, atom indices)."""
    names, resnames = _attr(u.atoms, 'names'), np.char.upper(_attr(u.atoms, 'resnames'))
    resix = u.atoms.resindices
    has = lambda *wanted: (np.bincount(resix[np.isin(names, wanted)], minlength=resix.max() + 1) > 0)[resix]
    protein, nucleic = has('N') & has('CA'), has("C4'", 'C4*')
    water, ion = np.isin(resnames, list(WATER)), np.isin(resnames, list(IONS))
    hydrogen = np.array([n.lstrip('0123456789')[:1] == 'H' for n in names])
    out = [('protein', protein), ('protein and nucleic acid', protein | nucleic),
           ('everything but water and ions', ~water & ~ion), ('everything but water', ~water),
           ('protein without hydrogens', protein & ~hydrogen),
           ('protein and nucleic acid without hydrogens', (protein | nucleic) & ~hydrogen),
           ('the backbone', protein & np.isin(names, ['N', 'CA', 'C', 'O'])),
           ('the Cα atoms', protein & (names == 'CA'))]
    return [(what, np.where(mask)[0]) for what, mask in out]


BONDED = (('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ("C4'", "C3'"), ("C1'", "C2'"), ('P', "O5'"),
          ('C4*', 'C3*'), ('C1*', 'C2*'), ('P', 'O5*'))


def _misfit(u) -> Optional[str]:
    """Do the coordinates fit the topology? Bonded atoms of a residue sit 1.2 to 1.6 Å apart. Returns None when
    they do, else the bond that does not. Several bonds are looked at, at both ends of the residue: programs
    order the atoms within a residue differently (GROMACS puts C and O last, OpenMM and the PDB before the side
    chain), so the first atoms alone can sit where they should while the rest do not."""
    names, resix = _attr(u.atoms, 'names'), u.atoms.resindices
    xyz = u.atoms.positions
    for one, two in BONDED:
        i, j = np.where(names == one)[0], np.where(names == two)[0]
        both = np.intersect1d(resix[i], resix[j])
        if len(both) >= 3:
            i, j = i[np.isin(resix[i], both)][:300], j[np.isin(resix[j], both)][:300]
            if len(i) == len(j):
                d = float(np.median(np.linalg.norm(xyz[i] - xyz[j], axis=1)))
                if not 1.0 < d < 1.9:
                    return f'the bonded atoms {one} and {two} are {d:.1f} Å apart in it'
    return None   # no bond to judge by


def _fits(u) -> bool:
    return _misfit(u) is None


def _reconcile(mda, u, traj: str, n_traj: int, notes: list):
    """Topology and trajectory with different atom counts: the trajectory was saved for part of the system
    (gmx trjconv/mdrun with an output group, cpptraj strip), or the topology was written for part of it."""
    import os
    n_top = len(u.atoms)
    tname = os.path.basename(traj)
    if n_traj < n_top:
        for what, idx in _groups(u):
            if len(idx) != n_traj:
                continue
            sub = mda.Merge(u.atoms[idx])
            sub.load_new(traj)
            if _fits(sub):
                notes.append(f'{tname} holds {n_traj} atoms, the topology {n_top}: the trajectory was taken to '
                             f'be {what} of the topology')
                return sub
        raise ValueError(f'the trajectory holds {n_traj} atoms and the topology {n_top}, and no usual part of '
                         'the topology (protein, protein and nucleic acid, everything but water) has that many. '
                         'Write the topology for the atoms the trajectory was saved for')
    # more atoms in the trajectory: the topology is taken to describe its first atoms (solute first, as MD
    # programs order a system)
    a = u.atoms
    chain_ids, icodes = _attr(a, 'chainIDs'), _attr(u.residues, 'icodes')
    pad = n_traj - n_top
    big = _empty(mda, list(a.names) + ['X'] * pad, list(a.resnames) + ['PAD'] * pad,
                 np.r_[a.resids, np.zeros(pad, int)], np.r_[a.segindices, np.full(pad, a.segindices.max() + 1)],
                 segids=np.r_[_attr(a, 'segids') if _attr(a, 'segids') is not None else np.full(n_top, ''),
                              np.full(pad, 'PAD')],
                 chain_ids=None if chain_ids is None else np.r_[chain_ids, np.full(pad, '')],
                 icodes=None if icodes is None else np.r_[icodes[a.resindices], np.full(pad, '')])
    big.load_new(traj)
    wrong = _misfit(big)
    if wrong:
        raise ValueError(f'the trajectory holds {n_traj} atoms and the topology {n_top}, and the topology does '
                         f'not describe the first atoms of the trajectory: {wrong}. Write the topology for the '
                         'same atoms as the trajectory')
    notes.append(f'{tname} holds {n_traj} atoms, the topology {n_top}: the topology was taken to describe the '
                 f'first {n_top} atoms of the trajectory (waters cannot be followed this way)')
    return big


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
        if f and os.path.isdir(f):
            raise ValueError(f'{f} is a folder, not a file')
        if f and not os.path.getsize(f):
            raise ValueError(f'{f} is empty')
    if _is_gromacs_top(top) and not traj:
        raise ValueError(f'{os.path.basename(top)} lists the atoms but holds no coordinates: give a trajectory '
                         'with it, or the .gro or .tpr instead')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # MDAnalysis is chatty about PDB records and future changes
        try:
            n_traj = _n_atoms(mda, traj) if traj else None
            u, notes = _topology(mda, top, n_traj)
            _original_numbers(mda, u, top, notes)
            if traj and n_traj != len(u.atoms):
                u = _reconcile(mda, u, traj, n_traj, notes)
            elif traj:
                u.load_new(traj)
                wrong = _misfit(u)
                if wrong:
                    raise ValueError(
                        f'the trajectory has as many atoms as the topology ({n_traj}) but does not belong to it: '
                        f'{wrong}. The two files list the atoms in a different order (written by different '
                        'programs?) or are different systems. Give the topology the trajectory was run with')
        except Exception as e:
            raise ValueError(f'Could not read {top}' + (f' with {traj}' if traj else '') +
                             f': {str(e).rstrip(".")}.') from e
    u._ppitango_notes = notes
    return u


def _polymer(u, renumber: Optional[dict] = None, extra=None):
    """Heavy atoms of protein and nucleic acid residues, with a chain label for every atom.

    Returns (atoms, index, table): `atoms` has one row per kept atom (chain, segid, resname, resnum, name,
    element), `index` is each row's position in the trajectory's atom order and `table` describes the chains.
    Residues count as polymer when they have a backbone (N and CA, or a sugar C4'), whatever they are called,
    so force-field names (HIE, CYX, ...) and modified residues stay and water, ions, lipids and ligands go,
    unless their atoms are listed in `extra` (a ligand, ion or sugar named in a selection).
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
    first = np.array([n[:1] for n in bare])
    # an element the atom name does not start with was guessed from a force-field atom type (a GROMACS .top
    # calls a histidine nitrogen NA, which reads as sodium): the name is the better guide
    given = elements if elements is not None else np.full(len(names), '')
    agrees = np.array([bool(e) and n.upper().startswith(e.upper()) for e, n in zip(given, bare)])
    elements = np.where(agrees | (np.char.upper(given) == 'D'), given, first)
    other = np.zeros(len(names), bool)
    if extra is not None and len(extra):
        other[extra] = True
        other &= ~(protein | nucleic)[resix]
        # outside a polymer the first letter is not enough: a lone atom is an ion (NA, ZN, MG), CL and BR halogens
        letters = np.array([''.join(ch for ch in n if ch.isalpha())[:2].upper() for n in names])
        lone = np.bincount(resix, minlength=n_res)[resix] == 1
        elements = np.where(other & ~agrees & (lone | np.isin(letters, ('CL', 'BR'))), letters, elements)
    hydrogen = np.isin(np.char.upper(elements), ('H', 'D'))
    altloc = _attr(ag, 'altLocs')
    keep = ((protein | nucleic)[resix] | other) & ~hydrogen
    if altloc is not None:
        keep &= np.isin(altloc, ('', 'A'))
    index = np.where(keep)[0]
    if not len(index):
        raise ValueError('No protein or nucleic acid found in the topology')

    ribose = has("O2'", 'O2*')[resix[index]]
    kind = np.where(protein[resix[index]], 'protein', np.where(nucleic[resix[index]], 'nucleic', 'other'))
    pairs = [canonical(r, n, rib) if k != 'other' else (r, n)
             for r, n, rib, k in zip(_attr(ag, 'resnames')[index], names[index], ribose, kind)]
    resname, name = np.array([q[0] for q in pairs]), np.array([q[1] for q in pairs])
    resnum = ag.resids[index]
    segid = _attr(ag, 'segids')
    segid = segid[index] if segid is not None else np.full(len(index), '')
    chain_id = _attr(ag, 'chainIDs')
    chain_id = chain_id[index] if chain_id is not None else np.full(len(index), '')
    icode = _attr(u.residues, 'icodes')
    icode = icode[resix[index]] if icode is not None else np.full(len(index), '')

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
            new = segid[s] != segid[s - 1] or resnum[s] < resnum[s - 1] or kind[s] != kind[s - 1] \
                or kind[s] == 'other'
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
                          'element': np.char.upper(elements[index]), 'kind': kind, 'icode': icode})
    for c, first in (renumber or {}).items():   # tpr and prmtop files count residues 1, 2, 3 ... through all chains
        if c not in set(chain):
            raise ValueError(f'--renumber: no chain {c}. Chains: {", ".join(dict.fromkeys(chain))}')
        rows = atoms.chain == c
        atoms.loc[rows, 'resnum'] += first - atoms.resnum[rows].iloc[0]
    rows = []
    for c, g in atoms.groupby('chain', sort=False):
        first, last = g.iloc[0], g.iloc[-1]
        n = g[['segid', 'resnum', 'icode']].drop_duplicates().shape[0]
        seg = f', segment {first.segid}' if first.segid and first.segid != 'SYSTEM' else ''
        what = {'protein': 'residues', 'nucleic': 'nucleotides'}.get(first.kind, 'other residues')
        rows.append(f'  {c}: {n} {what}, {first.resname}{first.resnum}{first.icode} to '
                    f'{last.resname}{last.resnum}{last.icode}{seg}')
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
    """Which atoms make up side a and side b. Returns (atoms, index, (a, b), log, sides): `sides` is what to
    call the two when they are not one chain each (several chains, a ligand, an ion), else ()."""
    log = list(getattr(u, '_ppitango_notes', []))
    chosen = []
    if select:
        if len(select) != 2:
            raise ValueError('Give exactly two selections, e.g. --select "resid 1-110" "resid 111-199"')
        for text in select:
            try:
                chosen.append(u.select_atoms(text).indices)
            except Exception as e:
                raise ValueError(f'Selection "{text}" not understood: {e}') from e
    # a selection that names no protein or nucleic acid at all is a ligand, a sugar, a lipid or ions: its heavy
    # atoms are followed as they are. Otherwise only the polymer in a selection counts, so "chainID A" does not
    # drag the waters and buffer molecules of chain A along.
    names, resix = _attr(u.atoms, 'names'), u.atoms.resindices
    backbone = np.isin(resix, resix[np.isin(names, ['CA', "C4'", 'C4*'])])
    extra = [found for found in chosen if len(found) and not backbone[found].any()]
    atoms, index, table = _polymer(u, renumber, extra=np.concatenate(extra) if extra else None)
    picked = None
    if select:
        pos = {j: i for i, j in enumerate(index)}
        polymer = (atoms.kind != 'other').values
        picked = [[pos[j] for j in found if j in pos and (polymer[pos[j]] or not backbone[found].any())]
                  for found in chosen]
        for rows, text in zip(picked, select):
            if not rows:
                raise ValueError(f'Selection "{text}" holds no atoms to follow (hydrogens are left out).\n{table}')
        if set(picked[0]) & set(picked[1]):
            raise ValueError('The two selections overlap')
        log.append(f'Side A = "{select[0]}", side B = "{select[1]}"')
    elif not chains and {'protein', 'nucleic'} <= set(atoms.kind):
        picked = [np.where(atoms.kind == k)[0].tolist() for k in ('protein', 'nucleic')]
        log += [table, 'Protein with nucleic acid: side A is all protein chains, side B all nucleic acid chains '
                       '(choose others with --chains or --select)']
    sides = ()
    if picked:
        rows = picked[0] + picked[1]
        atoms, index = atoms.iloc[rows].reset_index(drop=True), index[rows]
        atoms['part'] = atoms.chain          # the chain each atom came from: a side may hold several
        # in the frames PandaPPI reads the two sides are chains A and B, whatever the file calls their chains
        atoms['chain'] = ['A'] * len(picked[0]) + ['B'] * len(picked[1])
        a, b = 'A', 'B'
        sides = []
        for c in (a, b):
            side = atoms[atoms.chain == c]
            parts = list(dict.fromkeys(side.part))
            if (side.kind == 'other').all():     # a ligand or ions: by residue name
                names = list(dict.fromkeys(side.resname))
                sides.append(' + '.join(names) if len(names) <= 3 else f'{len(names)} kinds of residues')
            else:
                sides.append(f'Chain {parts[0]}' if len(parts) == 1 else
                             f"Chains {'+'.join(parts)}" if len(parts) <= 4 else f'{len(parts)} chains')
        sides = tuple(sides)
        log.append(f'Side A is {sides[0]} of the file, side B is {sides[1]}: the tables call them A and B')
    else:
        present = list(dict.fromkeys(atoms.chain))
        if chains:
            if len(chains) != 2 or chains[0] == chains[1]:
                raise ValueError('Give two different chains, e.g. --chains A D')
            a, b = chains
            missing = [c for c in (a, b) if c not in present]
            if missing:
                raise ValueError(f'Chain {" and ".join(missing)} not found.\n{table}\n'
                                 'Pick two of these with --chains, or describe the two sides with --select.')
        else:
            a, b = _closest_pair(atoms, u.atoms.positions[index].astype(float), u.trajectory.ts.dimensions)
            log.append(table)
            log.append(f'Following chains {a} and {b} (the pair with most contacts); choose others with --chains')
        if 'lettered' in table and chains:
            log.append(table)
        rows = np.r_[np.where(atoms.chain == a)[0], np.where(atoms.chain == b)[0]]
        atoms, index = atoms.iloc[rows].reset_index(drop=True), index[rows]
    atoms['unit'] = pd.factorize(atoms.part if 'part' in atoms else atoms.chain)[0]   # chains, for unwrapping
    atoms = _number(atoms)
    for c in (a, b):
        side = (atoms.chain == c).values
        if (atoms.resid[side].str.contains(':')).any():
            log.append(f'Side {c} holds several chains: its residues are named chain:number '
                       f'({atoms.resid[side].iloc[0]} ... {atoms.resid[side].iloc[-1]})')
    return atoms.drop(columns='segid'), index, (a, b), log, sides


def _whole(xyz: np.ndarray, unit: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Undo periodic wrapping. `unit` numbers the chains, side a's first: each chain is made whole, then put at
    the periodic image closest to the chains before it."""
    inv = np.linalg.inv(box)
    placed = None
    for k in dict.fromkeys(unit):
        sel = unit == k
        part = xyz[sel]
        # consecutive atoms of a chain are a few Å apart, so any jump of about a box length is a wrap
        part[1:] -= np.cumsum(np.round(np.diff(part, axis=0) @ inv), axis=0) @ box
        if placed is not None:
            ref, own = placed[::4], part[::4]
            base = np.round((own.mean(0) - ref.mean(0)) @ inv)
            best, shift = np.inf, base
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        s = base + (dx, dy, dz)
                        d = ((ref[:, None] - (own - s @ box)[None]) ** 2).sum(-1).min()
                        if d < best:
                            best, shift = d, s
            part -= shift @ box
        xyz[sel] = part
        placed = part if placed is None else np.concatenate([placed, part])
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

    Returns (atoms, coords, (a, b), times, log, waters, sides); `times` is in ps when the trajectory stores them,
    else None; `waters` is a list of oxygen positions per frame when asked for, else None.
    """
    if traj is None and not (select or renumber or waters) and \
            (_plain_pdb(top) or (chains and top.lower().endswith(('.pdb', '.ent')))):
        a, b = pick_chains(top, chains)
        atoms, coords, _ = read_pdb(top, (a, b), all_models=True)
        keep = np.arange(len(coords))[::stride]
        if max_frames and len(keep) > max_frames:
            keep = keep[np.linspace(0, len(keep) - 1, max_frames).round().astype(int)]
        order = np.r_[np.where(atoms.chain == a)[0], np.where(atoms.chain == b)[0]]
        return atoms.iloc[order].reset_index(drop=True), coords[keep][:, order], (a, b), None, [], None, ()
    import warnings
    u = _universe(top, traj)
    atoms, index, (a, b), log, sides = _sides(u, chains, select, renumber)
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
    b_sel, unit = (atoms.chain == b).values, atoms.unit.values
    frames, times, unwrapped = [], [], False
    for i in keep:
        ts = u.trajectory[i]
        xyz = ts.positions[index].astype(float)
        box = None
        if ts.dimensions is not None and ts.dimensions[:3].min() > 0:
            before = xyz.copy()
            box = ts.triclinic_dimensions.astype(float)
            xyz = _whole(xyz, unit, box)
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
    timed = bool(set(np.atleast_1d(u.trajectory.format)) & set(TIMED)) and len(times) > 1 and (np.diff(times) > 0).all()
    if waters:
        log.append(f'Interface waters: {np.mean([len(w) for w in found]):.1f} per frame on average sit within '
                   '3.5 Å of both sides')
    return atoms, np.array(frames), (a, b), times if timed else None, log, found, sides


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
    atoms, coords, (a, b), times, log, wat, sides = _read_trajectory(top, traj, chains, stride, max_frames,
                                                                     select, renumber, waters)
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
    return Ensemble(atoms, coords, (a, b), x=x, xlabel=xlabel, kind='md', waters=wat, log=log, sides=sides,
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
        atoms, coords, pair, _, log_k, wat_k, sides = _read_trajectory(
            w if is_pdb else top, None if is_pdb else w, chains, 1, frames_per_window, select, renumber, waters)
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
                    sides=sides,
                    title=f'{_name(top)} · free energy perturbation' + (f' · {label}' if label else ''),
                    curve=curve, curve_label='ΔΔG (kcal/mol)' if table is not None else '',
                    notes=[f'{len(windows)} λ windows, {len(coords)} frames',
                           'contacts are geometric: decoupled atoms still count'])
