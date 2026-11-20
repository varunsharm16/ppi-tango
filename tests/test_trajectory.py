"""Trajectory input the way MD programs write it: no chain IDs, wrapped across the periodic box, force-field
atom and residue names. Every variant must give the contacts the plain PDB gives."""
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
