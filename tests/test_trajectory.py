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


def test_side_of_two_chains_is_renumbered():
    ens = ensembles.md('examples/data/1BJ1.pdb', select=['chainID H L', 'chainID W'])
    assert any('shifted' in line for line in ens.log)
    side = ens.atoms[ens.atoms.chain == 'A']
    assert side.resnum.max() > 1000
