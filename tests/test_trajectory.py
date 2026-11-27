"""Trajectory input the way MD programs write it: no chain IDs, wrapped across the periodic box, force-field
atom and residue names. Every variant must give the contacts the plain PDB gives."""
import re

import numpy as np
import pytest

from ppitango import analysis, ensembles

mda = pytest.importorskip('MDAnalysis')
TOP = 'examples/simulations/output/md_1brs/top.pdb'
TRAJ = 'examples/simulations/output/md_1brs/traj.dcd'
FRAMES = 4


def _contacts(ens):
    """Every contact of every frame, in a fixed order (PandaPPI's own order depends on the chain letters)."""
    df = analysis.analyze(ens, workers=1, progress=False).contacts.fillna('')
    df = df.sort_values(['frame', 'type', 'residue_a', 'residue_b', 'atom_a', 'atom_b']).reset_index(drop=True)
    return df[['frame', 'type', 'residue_a', 'residue_b', 'atom_a', 'atom_b']], df.distance.values


@pytest.fixture(scope='module')
def reference():
    return _contacts(ensembles.md(TOP, TRAJ, chains=['A', 'D'], max_frames=FRAMES, dt=0.1))


@pytest.fixture(scope='module')
def wrapped(tmp_path_factory):
    """The same run as an MD program leaves it: a .gro without chain IDs and a trajectory with every atom
    wrapped into a box whose corner sits in the middle of the complex, so both chains are cut. (DCD rather than
    xtc: xtc rounds to 0.01 Å, enough to flip a contact sitting exactly on a cutoff.)"""
    import warnings
    out = tmp_path_factory.mktemp('gmx')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(TOP, TRAJ)
        keep = np.linspace(0, len(u.trajectory) - 1, FRAMES).round().astype(int)
        u.trajectory[0]
        u.atoms.write(str(out / 'top.gro'))
        with mda.Writer(str(out / 'traj.dcd'), u.atoms.n_atoms) as w:
            for ts in u.trajectory[keep]:
                box = ts.dimensions[:3]
                shifted = ts.positions - u.atoms.center_of_geometry()   # box corner at the interface
                ts.positions = shifted - np.floor(shifted / box) * box
                w.write(u.atoms)
    return str(out / 'top.gro'), str(out / 'traj.dcd')


def test_no_chain_ids_and_wrapped_box(reference, wrapped):
    ens = ensembles.md(*wrapped)                       # chains found by backbone breaks, pair picked by contact
    assert ens.chains == ('A', 'B')
    assert any('made whole' in line for line in ens.log)
    pairs, dist = _contacts(ens)
    ref_pairs, ref_dist = reference
    assert pairs.equals(ref_pairs)
    assert np.abs(dist - ref_dist).max() < 0.01


def test_select_two_sides(reference, wrapped):
    ens = ensembles.md(*wrapped, select=['resid 3-110 and index 0-1727', 'index 1728-9999'])
    assert _contacts(ens)[0].equals(reference[0])


def test_force_field_names(tmp_path, reference):
    """CHARMM spellings (HSD, ILE CD, OT1/OT2) and no chain column: same residues, same contacts."""
    lines = []
    for line in open(TOP):
        if line.startswith('ATOM'):
            name, res = line[12:16], line[17:20]
            res = 'HSD' if res == 'HIS' else res
            name = ' CD ' if (res, name) == ('ILE', ' CD1') else {' OXT': ' OT2'}.get(name, name)
            line = line[:12] + name + line[16] + res + '  ' + line[22:]
        lines.append(line)
    top = tmp_path / 'charmm.pdb'
    top.write_text(''.join(lines))
    ens = ensembles.md(str(top), TRAJ, max_frames=FRAMES, dt=0.1)
    assert set(ens.atoms.resname) >= {'HIS'} and 'HSD' not in set(ens.atoms.resname)
    assert _contacts(ens)[0].equals(reference[0])


def test_wrong_chain_lists_the_chains(wrapped):
    with pytest.raises(ValueError, match='Chain D not found(.|\n)*A: 108 residues'):
        ensembles.md(*wrapped, chains=['A', 'D'])


def test_side_of_two_chains_keeps_its_numbers():
    """Heavy and light chain on one side, their numbers repeating: the frames get running numbers, the residues
    keep their own, with the chain in front."""
    ens = ensembles.md('examples/data/1BJ1.pdb', select=['chainID H L', 'chainID W'])
    side = ens.atoms[ens.atoms.chain == 'A']
    assert side.resid.str.startswith(('H:', 'L:')).all() and {'H:1', 'H:100', 'L:1'} <= set(side.resid)
    assert side.resnum.max() == side.resid.nunique() <= 9999
    assert not ens.atoms[ens.atoms.chain == 'B'].resid.str.contains(':').any()


def test_renumber_restores_original_numbers(wrapped):
    """tpr and prmtop files count residues straight through; --renumber puts a chain back on its own numbers."""
    ens = ensembles.md(*wrapped, renumber={'A': 1, 'B': 501})
    assert ens.atoms.resnum[ens.atoms.chain == 'A'].min() == 1
    assert ens.atoms.resnum[ens.atoms.chain == 'B'].min() == 501
    with pytest.raises(ValueError, match='no chain Q'):
        ensembles.md(*wrapped, renumber={'Q': 1})


def test_force_field_nucleotide_names():
    assert ensembles.canonical('DA5', "O1P") == ('DA', 'OP1')
    assert ensembles.canonical('THY', "C5M") == ('DT', 'C7')
    assert ensembles.canonical('URA', "O2*", ribose=True) == ('U', "O2'")
    assert ensembles.canonical('RG3', 'N1') == ('G', 'N1')
    assert ensembles.canonical('HSD', 'OT1') == ('HIS', 'O')
    assert ensembles.canonical('ILE', 'CD') == ('ILE', 'CD1')
    assert ensembles.canonical('DT', 'O2') == ('DT', 'O2')      # a base oxygen, not a C-terminus


def test_phosphate_salt_bridge():
    pdb = ('ATOM      1  NZ  LYS A  10       0.000   0.000   0.000  1.00  0.00           N\n'
           'ATOM      2  OP1  DG B   5       3.000   0.000   0.000  1.00  0.00           O\n'
           'ATOM      3  OP2  DG B   6       9.000   0.000   0.000  1.00  0.00           O\n')
    assert analysis._extra_contacts(pdb, 'A', 'B') == [('Salt bridge / ionic', 'LYS 10', 'NZ', 'DG 5', 'OP1', 3.0)]
    assert analysis._extra_contacts(pdb, 'B', 'A')[0][1:3] == ('DG 5', 'OP1')   # side a always first


def test_interface_water_is_followed(tmp_path):
    """One water put between the two chains is picked up, moved with the frame and reported by PandaPPI."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(TOP)
    polar = u.select_atoms('(name N* or name O*) and not name H*')
    a, b = polar.select_atoms('chainID A'), polar.select_atoms('chainID D')
    d = np.linalg.norm(a.positions[:, None] - b.positions[None], axis=-1)
    i, j = np.argwhere((d > 5.0) & (d < 5.6))[0]
    w = (a.positions[i] + b.positions[j]) / 2
    lines = [l for l in open(TOP) if l.startswith(('CRYST1', 'ATOM'))]
    water = [f'HETATM{9000 + k:5d} {n:<4} HOH W   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {n[0]}\n'
             for k, (n, (x, y, z)) in enumerate([(' O', w), (' H1', w + 0.9), (' H2', w - 0.9)])]
    far = [f'HETATM 9900  O   HOH W   2    {w[0] + 30:8.3f}{w[1]:8.3f}{w[2]:8.3f}  1.00  0.00           O\n']
    top = tmp_path / 'wet.pdb'
    top.write_text(''.join(['MODEL        1\n'] + lines + water + far + ['ENDMDL\n'] +
                           ['MODEL        2\n'] + lines + water + far + ['ENDMDL\n']))
    ens = ensembles.md(str(top), chains=['A', 'D'], waters=True)
    assert [len(x) for x in ens.waters] == [1, 1]
    types = set(analysis.analyze(ens, workers=1, progress=False).contacts.type)
    assert 'Water-mediated' in types
    with pytest.raises(ValueError, match='no water'):
        ensembles.md(TOP, TRAJ, chains=['A', 'D'], waters=True, max_frames=2)


def test_gromacs_top(tmp_path, reference, wrapped):
    """A GROMACS .top as the topology: one molecule type per chain, and force-field atom types that look like
    other elements (the histidine nitrogen type NA is not sodium)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(wrapped[0])
    cut = int(np.where(np.diff(u.atoms.resids) < 0)[0][0]) + 1            # first atom of the second chain
    text = '[ defaults ]\n1 2 yes 0.5 0.8333\n\n'
    for mol, atoms in (('Protein_chain_A', u.atoms[:cut]), ('Protein_chain_D', u.atoms[cut:])):
        text += f'[ moleculetype ]\n{mol} 3\n\n[ atoms ]\n'
        for i, a in enumerate(atoms, 1):
            kind = 'NA' if a.name in ('ND1', 'NE2') else a.name[0]
            text += f'{i} {kind} {a.resid} {a.resname} {a.name} {i} 0.0 12.0\n'
        text += '\n'
    text += '[ system ]\ntest\n\n[ molecules ]\nProtein_chain_A 1\nProtein_chain_D 1\n'
    top = tmp_path / 'topol.top'
    top.write_text(text)
    ens = ensembles.md(str(top), wrapped[1])            # each molecule keeps its own residue numbers
    assert ens.chains == ('A', 'B')
    assert _contacts(ens)[0].equals(reference[0])
    with pytest.raises(ValueError, match='no coordinates'):
        ensembles.md(str(top), None)


def test_trajectory_saved_for_part_of_the_system(tmp_path, reference):
    """A trajectory written for the heavy atoms only, against the topology of everything: the part is recognised."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(TOP, TRAJ)
        heavy = u.select_atoms('not (name H* or name [123]H*)')
        with mda.Writer(str(tmp_path / 'heavy.dcd'), heavy.n_atoms) as w:
            for ts in u.trajectory[np.linspace(0, len(u.trajectory) - 1, FRAMES).round().astype(int)]:
                w.write(heavy)
    ens = ensembles.md(TOP, str(tmp_path / 'heavy.dcd'), chains=['A', 'D'])
    assert any('without hydrogens' in line for line in ens.log)
    assert _contacts(ens)[0].equals(reference[0])


def test_topology_for_the_first_atoms_only(tmp_path, reference):
    """The reverse: a solute-only topology with a trajectory that goes on (solvent after the solute)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        u = mda.Universe(TOP, TRAJ)
        n = u.atoms.n_atoms
        with mda.Writer(str(tmp_path / 'more.dcd'), n + 300) as w:
            big = mda.Universe.empty(n + 300, trajectory=True)
            for ts in u.trajectory[np.linspace(0, len(u.trajectory) - 1, FRAMES).round().astype(int)]:
                big.atoms.positions = np.r_[ts.positions, ts.positions[:300] + 40]
                big.dimensions = ts.dimensions
                w.write(big.atoms)
    ens = ensembles.md(TOP, str(tmp_path / 'more.dcd'), chains=['A', 'D'])
    assert any('first' in line for line in ens.log)
    assert _contacts(ens)[0].equals(reference[0])
    with pytest.raises(ValueError, match='no usual part'):       # a count nothing explains is still an error
        with mda.Writer(str(tmp_path / 'odd.dcd'), 777) as w:
            w.write(u.atoms[:777])
        ensembles.md(TOP, str(tmp_path / 'odd.dcd'))


def test_insertion_codes(tmp_path):
    """Residues 38 and 38A are two residues, and keep those names."""
    lines = []
    for line in open('examples/data/1BRS.pdb'):
        if line.startswith('ATOM') and line[21] == 'D' and int(line[22:26]) >= 39:
            num = int(line[22:26]) - 1
            line = line[:22] + f'{num:4d}' + ('A' if num == 38 else ' ') + line[27:]
        lines.append(line)
    pdb = tmp_path / 'icode.pdb'
    pdb.write_text(''.join(lines))
    for kw in (dict(chains=['A', 'D']), dict(select=['chainID A', 'chainID D'])):    # both readers
        ens = ensembles.nma(str(pdb), frames=2, **kw)
        side = ens.atoms[ens.atoms.chain == ens.chains[1]]
        assert {'38', '38A', '39'} <= set(side.resid)
        assert side.groupby('resid').resnum.nunique().max() == 1 and side.resnum.nunique() == side.resid.nunique()
    scene = analysis.analyze(ens, workers=1, progress=False).scene()
    assert any(c['residue_b'].endswith(' 38A') for c in scene.contacts)


def test_metal_ion_as_a_side():
    """A selection without protein in it (here the zinc of ACE2) is followed as it is."""
    ens = ensembles.nma('examples/data/6M0J.pdb', select=['chainID A', 'resname ZN'], frames=2)
    assert ens.atoms[ens.atoms.kind == 'other'].element.tolist() == ['ZN']
    assert ens.sides == ('Chain A', 'ZN') and ens.pair == 'chain A and ZN'   # named as the file names them
    contacts = analysis.analyze(ens, workers=1, progress=False).contacts
    assert (contacts.type == 'Metal coordination').sum() >= 3


def test_base_stacking_and_ligand_contacts():
    hexagon = [(1.4 * np.cos(k * np.pi / 3), 1.4 * np.sin(k * np.pi / 3)) for k in range(6)]
    row = lambda i, name, res, chain, num, x, y, z, el: (
        f'ATOM  {i:5d} {name:<4} {res:>3} {chain}{num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}\n')
    pdb = ''.join(row(i, n, 'PHE', 'A', 7, x, y, 0, 'C')
                  for i, (n, (x, y)) in enumerate(zip(('CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2'), hexagon)))
    pdb += ''.join(row(9 + i, n, 'DC', 'B', 3, x, y, 3.6, n[0])
                   for i, (n, (x, y)) in enumerate(zip(('N1', 'C2', 'N3', 'C4', 'C5', 'C6'), hexagon)))
    assert [r[:2] for r in analysis._extra_contacts(pdb, 'A', 'B')] == [('Pi-pi stacking', 'PHE 7')]
    pdb = row(1, 'CD1', 'LEU', 'A', 5, 0, 0, 0, 'C') + row(2, 'C7', 'LIG', 'B', 1, 4, 0, 0, 'C') + \
        row(3, 'C8', 'LIG', 'B', 1, 4, 1.4, 0, 'C') + row(4, 'O1', 'LIG', 'B', 1, 4, 2.6, 0, 'O')
    assert analysis._extra_contacts(pdb, 'A', 'B') == [('Hydrophobic', 'LEU 5', 'CD1', 'LIG 1', 'C7', 4.0)]


def test_prmtop_residue_numbers(tmp_path):
    """An Amber prmtop that carries the PDB numbers and chains (parmed addPDB) is read with them."""
    top = tmp_path / 'x.prmtop'
    top.write_text('%VERSION\n%FLAG RESIDUE_NUMBER\n%FORMAT(20I4)\n   3   4   1\n'
                   '%FLAG RESIDUE_CHAINID\n%FORMAT(20a4)\nA   A   D   \n')
    u = ensembles._empty(mda, ['N', 'CA', 'N', 'CA', 'N', 'CA'], ['ALA'] * 6, [1, 1, 2, 2, 3, 3], [0] * 6)
    notes = []
    ensembles._original_numbers(mda, u, str(top), notes)
    assert u.residues.resids.tolist() == [3, 4, 1] and u.atoms.chainIDs.tolist() == list('AAAADD') and notes


def test_trajectory_of_another_setup_is_refused(tmp_path):
    """The same protein written by another program has its atoms in another order within each residue (here C
    and O moved behind the side chain, as GROMACS has them): same atom count, wrong atoms."""
    lines = [l for l in open(TOP) if l.startswith(('ATOM', 'HETATM'))]
    out, residue = [], []
    for line in lines + ['']:
        if residue and line[17:27] != residue[0][17:27]:
            last = [l for l in residue if l[12:16].strip() in ('C', 'O')]
            out += [l for l in residue if l not in last] + last
            residue = []
        residue.append(line)
    assert len(out) == len(lines) and out != lines
    other = tmp_path / 'other.pdb'
    other.write_text(''.join(out) + 'END\n')
    with pytest.raises(ValueError, match='as many atoms as the topology.*does not belong to it'):
        ensembles.md(str(other), TRAJ, chains=['A', 'D'], max_frames=FRAMES)


def _renamed(path, change, only=None, at=None):
    """TOP with the atom names of every residue (of the kinds in `only`, or the one residue `at`, chain and
    number as columns 22 to 26 have them) handed round by `change`: the coordinates of the trajectory then land
    on other atoms of the same residue."""
    out, residue = [], []
    for line in [l for l in open(TOP) if l.startswith(('ATOM', 'HETATM'))] + ['']:
        if residue and line[17:27] != residue[0][17:27]:
            names = [(l[12:16], l[76:78]) for l in residue]
            hit = (only is None or residue[0][17:20] in only) and (at is None or residue[0][21:26] == at)
            new = change(names) if hit else names
            out += [l[:12] + a + l[16:76] + e + l[78:] for l, (a, e) in zip(residue, new)]
            residue = []
        residue.append(line)
    path.write_text(''.join(out) + 'END\n')
    return str(path)


def _swap(one, two):
    def change(names):
        bare = [n.strip() for n, _ in names]
        i, j = bare.index(one), bare.index(two)
        names = list(names)
        names[i], names[j] = names[j], names[i]
        return names
    return change


def _reverse_side_chain(names):
    """Heavy atoms after CB in reverse order, backbone untouched."""
    bare = [n.strip() for n, _ in names]
    if 'CB' not in bare:
        return names
    side = [k for k in range(bare.index('CB') + 1, len(bare))
            if not bare[k].startswith('H') and bare[k] not in ('C', 'O', 'OXT')]
    names = list(names)
    for k, value in zip(side, [names[k] for k in reversed(side)]):
        names[k] = value
    return names


@pytest.mark.parametrize('change, only, at, frames, said', [
    (_reverse_side_chain, None, None, FRAMES, 'in frame 1 the bonded atoms CB and CG of .* more bonds are off'),
    (_swap('NE', 'NH1'), ('ARG',), None, FRAMES, 'in frame 1 the bonded atoms .* of ARG'),
    # 0.1 Å between the two bonds: with four frames all residues of the kind together tell
    (_swap('OD1', 'ND2'), ('ASN',), None, FRAMES, r'over the 4 frames read from 1 to \d+, in the \d+ ASN, .*OD1 and ND2 have'),
    (_swap('OE1', 'NE2'), ('GLN',), None, FRAMES, r'in the \d+ GLN, .*OE1 and NE2 have changed places'),
    (_swap('OG1', 'CG2'), ('THR',), None, FRAMES, 'OG1 and CG2 have changed places'),
    # one residue only, at the interface and far from it: every residue is judged on its own
    (_reverse_side_chain, None, 'A  62', FRAMES, 'CB and CG of LYS 62 of chain A are'),
    (_reverse_side_chain, None, 'D  54', FRAMES, 'of ARG 54 of chain D are'),
    (_swap('CD1', 'CD2'), None, 'A  35', FRAMES, 'of TRP 35 of chain A are'),
    (_swap('CB', 'OG'), None, 'A  38', FRAMES, 'CA and CB of SER 38 of chain A are'),
    (_swap('C', 'O'), None, 'A  59', FRAMES, 'CA and C of ARG 59 of chain A are'),
    (_swap('N', 'CA'), None, 'A  27', FRAMES, 'of LYS 27 of chain A are'),
    (_swap('OD1', 'ND2'), None, 'A  58', 20, r'over the 20 frames read from 1 to \d+, in ASN 58 of chain A, .*OD1 and ND2'),
    (_swap('OE1', 'NE2'), None, 'D  72', 20, 'in GLN 72 of chain D, .*OE1 and NE2'),
    (_swap('OG1', 'CG2'), None, 'A  26', 20, 'in THR 26 of chain A, .*OG1 and CG2'),
])
def test_atoms_in_another_order_are_refused(tmp_path, change, only, at, frames, said):
    """The first frame's backbone as a whole looks right: only the bonds of each residue tell."""
    other = _renamed(tmp_path / 'other.pdb', change, only, at)
    with pytest.raises(ValueError, match=f'does not belong to .*: .*{said}'):
        ensembles.md(other, TRAJ, chains=['A', 'D'], max_frames=frames)


def _rewritten(path, change, n=24):
    """The first `n` frames of TRAJ with `change(positions, frame)` applied to each."""
    u = mda.Universe(TOP, TRAJ)
    with mda.Writer(str(path), len(u.atoms)) as w:
        for ts in u.trajectory[:n]:
            ts.positions = change(u, ts)
            w.write(u.atoms)
    return str(path)


def test_trajectory_that_turns_wrong_later_is_refused(tmp_path):
    """Two runs with differently ordered topologies joined into one file: right at first, wrong from frame 13
    on. Every frame that is read is measured, not a sample."""
    def change(u, ts):
        xyz = ts.positions.copy()
        if ts.frame >= 12:
            side = u.select_atoms('chainID A and resid 62 and name CG CD CE NZ').indices
            xyz[side] = xyz[side[::-1]]
        return xyz
    joined = _rewritten(tmp_path / 'joined.dcd', change)
    with pytest.raises(ValueError, match='in frame 13 the bonded atoms .* of LYS 62 .*The frames before fit'):
        ensembles.md(TOP, joined, chains=['A', 'D'], max_frames=None)

    def backbone(u, ts):     # the backbone is measured in every frame too, not in the first alone
        xyz = ts.positions.copy()
        if ts.frame >= 1:
            c, o = u.select_atoms('protein and name C').indices, u.select_atoms('protein and name O').indices
            xyz[c], xyz[o] = ts.positions[o], ts.positions[c]
        return xyz
    joined = _rewritten(tmp_path / 'backbone.dcd', backbone, 6)
    with pytest.raises(ValueError, match='in frame 2 the bonded atoms CA and C of .* more bonds are off'):
        ensembles.md(TOP, joined, chains=['A', 'D'], max_frames=None)


def test_noisy_coordinates_are_not_taken_for_a_wrong_trajectory(tmp_path):
    """A correct run stored coarsely (0.1 Å, xtc precision 100) or with 0.05 Å of noise on every coordinate
    spreads the bond lengths: the Asn, Gln and Thr rule asks for more in step, and no bond leaves its 0.4 Å."""
    rng = np.random.default_rng(1)
    for name, change in (('coarse', lambda u, ts: np.round(ts.positions, 1)),
                         ('noisy', lambda u, ts: ts.positions + rng.normal(0, 0.05, ts.positions.shape))):
        path = _rewritten(tmp_path / f'{name}.dcd', change, 200)
        assert ensembles.md(TOP, path, chains=['A', 'D'], max_frames=None).n_frames == 200
    path = _rewritten(tmp_path / 'spot.dcd', lambda u, ts: np.round(ts.positions), 2)
    with pytest.raises(ValueError, match='spot.dcd: in frame 1 .*stored with too few digits to tell atoms apart') as e:
        ensembles.md(TOP, path, chains=['A', 'D'])
    assert 'does not belong' not in str(e.value)


def test_short_reads_of_a_correct_run_are_never_refused():
    """One residue gets no verdict on the Asn, Gln and Thr rule from a handful of frames: a mean over five
    is too easily on the wrong side when the frames are close in time or coarsely stored."""
    u = ensembles._universe(TOP, TRAJ)
    rng = np.random.default_rng(0)
    for _ in range(150):
        n, stride = int(rng.integers(1, 30)), int(rng.integers(1, 4))
        start = int(rng.integers(0, len(u.trajectory) - n * stride))
        check = ensembles._Bonds(u)
        for frame in range(start, start + n * stride, stride):
            ts = u.trajectory[frame]
            ts.positions = np.round(ts.positions, 1)
            assert check.frame(ts) is None
        assert check.close() is None
    tally = ensembles._Tally(2)       # two frames that happen to agree: no spread, and still no verdict
    for frame in (1, 2, 3):
        tally.add(np.array([0.06, -0.09]), frame)
    assert not len(tally.above(0.045, 6.0, 20)[1]) and len(tally.above(0.045, 6.0, 3)[1]) == 1


def test_whole_residues_that_changed_places_are_refused(tmp_path):
    """Two lysines with each other's coordinates: every bond within a residue is right, the peptide bonds
    to their neighbours are not."""
    def change(u, ts):
        xyz = ts.positions.copy()
        a, b = (u.select_atoms(sel).indices for sel in ('chainID A and resid 27', 'chainID D and resid 21'))
        xyz[a], xyz[b] = ts.positions[b], ts.positions[a]
        return xyz
    path = _rewritten(tmp_path / 'exchanged.dcd', change, 3)
    said = 'the backbone N of LYS 27 of chain A is bonded to the C of LEU 20 of chain D, which is not the residue'
    with pytest.raises(ValueError, match=f'in frame 1 {said} before it .*whole residues have changed places'):
        ensembles.md(TOP, path, chains=['A', 'D'])


def test_amide_swap_in_part_of_a_run_is_refused(tmp_path):
    """An Asn at the interface with its oxygen and nitrogen exchanged in the second half of the frames: the
    mean over the whole run is near zero, the mean over 50 frames in a row is not."""
    def change(u, ts):
        xyz = ts.positions.copy()
        if ts.frame >= 100:
            o, n = (u.select_atoms(f'chainID D and resid 33 and name {x}').indices for x in ('OD1', 'ND2'))
            xyz[o], xyz[n] = ts.positions[n], ts.positions[o]
        return xyz
    path = _rewritten(tmp_path / 'half.dcd', change, 200)
    said = 'over the 50 frames read from 101 to 150, in ASN 33 of chain D, .*The frames before fit the topology'
    with pytest.raises(ValueError, match=said):
        ensembles.md(TOP, path, chains=['A', 'D'], max_frames=None)


def test_residues_cut_by_the_box_are_not_taken_for_a_wrong_trajectory(tmp_path, reference):
    """A raw trajectory has every atom put into the box on its own. With a box face through each of the three
    His (all there are), their bonds reach across the box in every frame: still the right trajectory."""
    def change(u, ts):
        box = ts.dimensions[:3]
        cut = [u.select_atoms(f'chainID {c} and resid {r} and name CG').positions[0][k]
               for k, (c, r) in enumerate((('A', 102), ('D', 17), ('A', 18)))]
        xyz = ts.positions - np.array(cut) + box
        return xyz - np.floor(xyz / box) * box
    u = mda.Universe(TOP, TRAJ)
    keep = np.linspace(0, len(u.trajectory) - 1, FRAMES).round().astype(int)
    with mda.Writer(str(tmp_path / 'raw.dcd'), len(u.atoms)) as w:
        for ts in u.trajectory[keep]:
            ts.positions = change(u, ts)
            w.write(u.atoms)
    ens = ensembles.md(TOP, str(tmp_path / 'raw.dcd'), chains=['A', 'D'], dt=0.1)
    assert any('made whole' in line for line in ens.log)
    keys, dist = _contacts(ens)       # the shift costs a float32 its last digits: distances are printed to 0.01 Å
    assert keys.equals(reference[0]) and np.abs(dist - reference[1]).max() < 0.02


def test_force_field_files_pass_the_side_chain_check():
    """GROMACS, Amber and CHARMM topologies name residues and atoms their own way (HIE, CYX, ILE CD)."""
    for top, traj in (('gromacs/md.gro', 'gromacs/md.xtc'), ('amber/complex.prmtop', 'amber/md.nc'),
                      ('charmm/complex.psf', 'charmm/md.dcd')):
        u = ensembles._universe(f'examples/formats/{top}', f'examples/formats/{traj}')
        check = ensembles._Bonds(u)
        assert len(check.i) > 1300 and len(check.pairs) == 3
        assert [check.frame(ts) for ts in u.trajectory] == [None] * len(u.trajectory) and check.close() is None


def test_missing_residues_do_not_split_a_chain(tmp_path):
    """A crystal structure without chain IDs: barstar lacks residues 64 and 65 in two of its three copies. The
    hole is one the two missing residues span, so it is one chain, and the pair followed is a real complex."""
    lines = [l[:21] + ' ' + l[22:] if l.startswith(('ATOM', 'HETATM', 'TER')) and len(l) > 22 else l
             for l in open('examples/data/1BRS.pdb')]
    for name, kept in (('ter.pdb', lines), ('no_ter.pdb', [l for l in lines if not l.startswith('TER')])):
        (tmp_path / name).write_text(''.join(kept))
        ens = ensembles.nma(str(tmp_path / name), frames=4)
        table = '\n'.join(ens.log)
        assert table.count(' residues, ') == 6 and 'D has no residues between THR63 and GLY66' in table
        assert {ens.atoms.resname[ens.atoms.chain == c].iloc[0] for c in ens.chains} == {'ALA', 'LYS'}
    # chain IDs on some atoms only: said so, and a chain is not cut where its ID stops
    some = [l[:21] + ' ' + l[22:] if l.startswith('ATOM') and l[21] == 'A' and int(l[22:26]) >= 50 else l
            for l in open('examples/data/1BRS.pdb')]
    (tmp_path / 'some.pdb').write_text(''.join(some))
    table = '\n'.join(ensembles.nma(str(tmp_path / 'some.pdb'), frames=4).log)
    assert 'some atoms have no chain ID' in table and 'A: 108 residues, VAL3 to ARG110' in table
    # a long missing stretch is one too, as long as its residues could span it (no limit on their number)
    long = [l for l in lines if not (l.startswith('ATOM') and 20 <= int(l[22:26]) <= 60
                                     and l[6:11] < ' 1000' and not l.startswith('TER'))]
    (tmp_path / 'long.pdb').write_text(''.join(long))
    table = '\n'.join(ensembles.nma(str(tmp_path / 'long.pdb'), frames=4).log)
    assert 'A has no residues between ' in table and table.count(' residues, ') == 6
    # one chain ID for everything, as MDAnalysis writes a selection without chains, is no chain ID
    same = [l[:21] + 'X' + l[22:72] + 'SYST' + l[76:] if l.startswith(('ATOM', 'HETATM')) else l for l in open(TOP)]
    (tmp_path / 'x.pdb').write_text(''.join(same))
    ens = ensembles.nma(str(tmp_path / 'x.pdb'), frames=4)
    assert ens.chains == ('A', 'B') and 'no chain IDs in the topology' in '\n'.join(ens.log)


def test_calpha_only_model_is_refused_with_the_reason(tmp_path):
    only = tmp_path / 'ca.pdb'
    only.write_text(''.join(l for l in open(TOP) if not l.startswith('ATOM') or l[12:16].strip() == 'CA'))
    with pytest.raises(ValueError, match='Cα-only and coarse-grained models are not read'):
        ensembles.nma(str(only), frames=4)
    with pytest.raises(ValueError, match='Cα-only and coarse-grained models are not read'):
        ensembles.nma(str(only), chains=['A', 'D'], frames=4)
    backbone = ('N', 'CA', 'C', 'O', 'CB')
    only.write_text(''.join(l for l in open(TOP) if not l.startswith('ATOM') or l[12:16].strip() in backbone))
    with pytest.raises(ValueError, match='165 of 165 residues have no side-chain atoms: this is a backbone model'):
        ensembles.nma(str(only), frames=4)
    some = [l for l in open(TOP) if not (l.startswith('ATOM') and l[21] == 'A' and 20 <= int(l[22:26]) < 30
                                         and l[12:16].strip() not in backbone and l[76:78].strip() != 'H')]
    only.write_text(''.join(some))
    log = ' '.join(ensembles.nma(str(only), frames=4).log)
    assert re.search(r'\d of 165 residues have no side-chain atoms in the file \(A:2\d, ', log)
    renamed = [l[:17] + 'ALA' + l[20:] for l in open(TOP) if l.startswith('ATOM') and l[12:16].strip() in backbone]
    only.write_text(''.join(renamed))
    assert 'Every amino acid is Gly, Ala or of an unknown kind' in ' '.join(ensembles.nma(str(only), frames=4).log)


def test_pdb_without_chain_ids_gets_its_chains_lettered(tmp_path):
    """`nma` and a multi-model PDB through `md` letter the chains as a topology with a trajectory does."""
    lines = [l[:21] + ' ' + l[22:] if l.startswith(('ATOM', 'TER')) and len(l) > 22 else l for l in open(TOP)]
    bare = tmp_path / 'bare.pdb'
    bare.write_text(''.join(lines))
    ens = ensembles.nma(str(bare), frames=4)
    assert ens.chains == ('A', 'B') and any('lettered' in line for line in ens.log)
    assert ensembles.nma(str(bare), chains=['A', 'B'], frames=4).coords.shape == ens.coords.shape


def test_netcdf_times_are_read():
    ens = ensembles.md('examples/formats/amber/complex.prmtop', 'examples/formats/amber/md.nc', max_frames=3)
    assert ens.xlabel == 'Time (ps)' and ens.x.tolist() == [0, 10, 20]
    assert ens.sides == ('Chain A', 'Chain D') and ens.pair == 'chains A and D'


def test_ring_contacts_end_on_the_ring_centre():
    """Cation-π, CH-π, π-π and amide-aromatic contacts name no atom on the ring's side and give the distance to
    its centre: the scene draws them there, so every drawn distance is the one in the table."""
    ens = ensembles.md(TOP, TRAJ, chains=['A', 'D'], max_frames=FRAMES)
    result = analysis.analyze(ens, workers=1, progress=False)
    scene = result.scene()
    table = result.contacts.sort_values('distance').drop_duplicates(['frame', 'type', 'residue_a', 'residue_b'])
    rings = table[(table.atom_a == '') | (table.atom_b == '')]
    assert set(rings.type) >= {'Cation-pi', 'CH-pi'} and scene.points.shape[1] > len(ens.atoms)
    index = {(c['type'], c['residue_a'], c['residue_b']): k for k, c in enumerate(scene.contacts)}
    for row in table.itertuples():
        k = index[(row.type, row.residue_a, row.residue_b)]
        assert scene.dist[k, row.frame] == pytest.approx(row.distance, abs=0.011)
    centres = np.unique(np.r_[scene.atom_a.ravel(), scene.atom_b.ravel()])
    assert (centres >= len(ens.atoms)).sum() == scene.points.shape[1] - len(ens.atoms)
