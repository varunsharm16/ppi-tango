"""End-to-end checks on a small synthetic two-chain complex (no downloads, a few seconds)."""
import os

import numpy as np
import pytest

import ppitango
from ppitango import analysis, ensembles, movie, viewer

# two short antiparallel strands, 4.8 Å apart, with side chains that can pair up
RESIDUES_A = ['LYS', 'SER', 'ARG', 'THR', 'ASN', 'TYR', 'GLN', 'ARG']
RESIDUES_B = ['GLU', 'THR', 'ASP', 'SER', 'GLN', 'PHE', 'ASN', 'GLU']


def _strand(chain, names, y, flip, serial):
    lines = []
    for i, name in enumerate(names):
        x = 3.8 * i
        side = 1 if (i % 2 == 0) != flip else -1
        atoms = [('N', x - 1.2, y, 0.3, 'N'), ('CA', x, y, 0.0, 'C'), ('C', x + 1.2, y, 0.3, 'C'),
                 ('O', x + 1.3, y + 0.9 * (-1 if flip else 1), 1.0, 'O'),
                 ('CB', x, y + 1.5 * (-1 if flip else 1), -0.5 * side, 'C')]
        tip = {'LYS': ('NZ', 'N'), 'ARG': ('NH1', 'N'), 'GLU': ('OE1', 'O'), 'ASP': ('OD1', 'O'),
               'SER': ('OG', 'O'), 'THR': ('OG1', 'O'), 'ASN': ('ND2', 'N'), 'GLN': ('NE2', 'N'),
               'TYR': ('OH', 'O'), 'PHE': ('CZ', 'C')}[name]
        atoms.append((tip[0], x, y + 2.6 * (-1 if flip else 1), 0.0, tip[1]))
        for an, ax, ay, az, el in atoms:
            serial += 1
            az += 0.9 * side + 0.06 * i * i      # pleat and bend the strand: a flat network has no stiffness
            lines.append(f'ATOM  {serial:5d}  {an:<3} {name} {chain}{i + 1:4d}    {ax:8.3f}{ay:8.3f}{az:8.3f}'
                         f'  1.00  0.00          {el:>2}')
    return lines, serial


@pytest.fixture(scope='module')
def pdb(tmp_path_factory):
    a, n = _strand('A', RESIDUES_A, 0.0, False, 0)
    b, _ = _strand('B', RESIDUES_B, 8.0, True, n)
    path = tmp_path_factory.mktemp('data') / 'toy.pdb'
    path.write_text('\n'.join(a + b) + '\nEND\n')
    return str(path)


@pytest.fixture(scope='module')
def result(pdb):
    ens = ensembles.nma(pdb, frames=6, rmsd=0.5)
    return analysis.analyze(ens, workers=2, progress=False)


def test_nma_oscillates_about_the_structure(pdb):
    ens = ensembles.nma(pdb, frames=8, rmsd=1.0)
    assert ens.chains == ('A', 'B') and ens.periodic
    assert ens.coords.shape[0] == 8
    start, _, _ = ensembles.read_pdb(pdb)
    assert np.allclose(ens.coords[0], ens.coords[4], atol=1e-6)          # phase 0 and 0.5: undisplaced
    ca = (ens.atoms.name == 'CA').values
    rmsd = np.sqrt(((ens.coords[2, ca] - ens.coords[0, ca]) ** 2).sum(-1).mean())
    assert rmsd == pytest.approx(1.0, abs=1e-6)                          # quarter phase: full swing
    assert len(start) == len(ens.atoms)


def test_contacts_tracked_through_frames(result):
    assert not result.contacts.empty
    assert set(result.contacts.frame) == set(range(6))
    scene = result.scene()
    assert scene.occ.shape == (len(scene.contacts), 6)
    assert scene.occ.min() >= 0 and scene.occ.max() <= 1
    assert scene.drawn and scene.left and scene.right
    # the distance of a contact is measured between its own two ends in every frame
    k = scene.drawn[0]
    xyz = scene.points
    d = np.linalg.norm(xyz[3, scene.atom_a[k, 3]] - xyz[3, scene.atom_b[k, 3]])
    assert scene.dist[k, 3] == pytest.approx(d)


def test_smoothing_turns_flicker_into_a_fade():
    flicker = np.array([[1, 0, 1, 0, 1, 1, 1, 0, 0, 0]], float)
    smooth = analysis._smooth(flicker, 5, periodic=False)[0]
    assert 0 < smooth[1] < 1 and smooth[-1] < smooth[5]
    assert np.allclose(analysis._smooth(flicker, 1, False), flicker)


def test_md_reads_multimodel_pdb(pdb, tmp_path):
    atoms, coords, _ = ensembles.read_pdb(pdb)
    path = tmp_path / 'traj.pdb'
    models = [f'MODEL {i + 1}\n' + ensembles._pdb_text(atoms, coords[0] + 0.1 * i).replace('END\n', 'ENDMDL\n')
              for i in range(4)]
    path.write_text(''.join(models))
    ens = ensembles.md(str(path), dt=0.5)
    assert ens.n_frames == 4 and ens.x[-1] == pytest.approx(1.5)
    assert np.allclose(ens.coords[3], ens.coords[0], atol=1e-3)          # rigid shifts are fitted away


def test_outputs(result, tmp_path):
    scene = result.scene()
    html = viewer.write(scene, str(tmp_path / 'v.html'))
    text = open(html, encoding='utf-8').read()
    assert '__DATA__' not in text and '$3Dmol' in text and '"contacts":[' in text
    mp4 = movie.render(scene, str(tmp_path / 'm.mp4'), seconds=0.4, fps=10)
    assert open(mp4, 'rb').read(12)[4:8] == b'ftyp'
    files = result.save(str(tmp_path / 'out'))
    assert all(open(f).readline().startswith('frame,x,') for f in files)


def test_cli_rejects_single_chain(pdb, tmp_path):
    lines = [l for l in open(pdb) if l[21:22] == 'A']
    one = tmp_path / 'one.pdb'
    one.write_text(''.join(lines))
    from ppitango.cli import main
    with pytest.raises(SystemExit) as e:
        main(['nma', str(one)])
    assert 'single chain' in str(e.value)


def test_cli_checks_options_before_working(pdb, tmp_path, capsys):
    from ppitango.cli import main
    for bad in (['--mode', '0'], ['--fps', '0'], ['--seconds', '0'], ['--frames', '1'], ['--workers', '0']):
        with pytest.raises(SystemExit) as e:
            main(['nma', pdb, *bad])
        assert e.value.code == 2 and 'must be' in capsys.readouterr().err
    for bad, said in ((['--chains', 'A', 'A'], 'two different chains'), (['--mode', '1000'], 'has 42 modes'),
                      (['-o', f'{pdb}/x'], 'cannot write to')):
        with pytest.raises(SystemExit) as e:
            main(['nma', pdb, *bad])
        assert said in str(e.value)
    for bad in (['--seconds', 'nan'], ['--rmsd', 'inf'], ['--distance-cutoff', 'nan']):
        with pytest.raises(SystemExit) as e:
            main(['nma', pdb, *bad])
        err = capsys.readouterr().err
        assert e.value.code == 2 and err.startswith('ppitango: ') and err.count('\n') == 1
    with pytest.raises(SystemExit) as e:
        main(['md', pdb, '--unit', 'parsecs'])
    assert e.value.code == 2 and 'invalid choice' in capsys.readouterr().err
    with pytest.raises(SystemExit) as e:
        main(['fep', pdb, pdb, '--dg', str(tmp_path / 'none.csv')])
    assert 'none.csv' in str(e.value)
    with pytest.raises(SystemExit) as e:
        main(['nma', pdb, '--renumber', 'A:x'])
    assert e.value.code == 2 and 'CHAIN:FIRST' in capsys.readouterr().err
    assert os.listdir(tmp_path) == []
    for text, said in (('window,dG\n0,abc\n1,def\n', "row 1 holds 'abc'"), ('window,dG\n0,0.5\n1,\n', 'row 2')):
        (tmp_path / 'dg.csv').write_text(text)
        with pytest.raises(SystemExit) as e:
            main(['fep', pdb, pdb, pdb, '--dg', str(tmp_path / 'dg.csv'), '-o', str(tmp_path / 'x')])
        assert 'dg.csv: the dG column needs a number' in str(e.value) and said in str(e.value)
    (tmp_path / 'dg.csv').write_bytes(b'A\xc5\x01\xff\n' * 20)
    with pytest.raises(SystemExit) as e:
        main(['fep', pdb, pdb, pdb, '--dg', str(tmp_path / 'dg.csv'), '-o', str(tmp_path / 'x')])
    assert 'dg.csv: not a text file' in str(e.value)


def test_cli_output_folder(pdb, tmp_path, capsys):
    """-o takes a prefix or a folder; a folder that cannot be written stops the run before the analysis."""
    from ppitango.cli import main
    name = os.path.splitext(os.path.basename(pdb))[0]
    main(['nma', pdb, '--frames', '4', '--no-movie', '--no-viewer', '-o', f'{tmp_path}/new/'])
    assert sorted(os.listdir(tmp_path / 'new')) == [f'{name}_nma_buried_surface.csv', f'{name}_nma_contacts.csv']
    main(['nma', pdb, '--frames', '4', '--no-movie', '--no-viewer', '-o', str(tmp_path / 'new')])
    assert len(os.listdir(tmp_path / 'new')) == 2 and os.listdir(tmp_path) == ['new']
    # an output file that cannot be replaced stops the run before the analysis too
    (tmp_path / 'new' / 'taken_viewer.html').mkdir()
    movie = tmp_path / 'new' / 'taken_contacts.mp4'
    movie.write_text('')
    movie.chmod(0o444)
    for options, said in ((['--no-movie'], 'taken_viewer.html: a folder of that name is in the way'),
                          (['--no-viewer'], 'taken_contacts.mp4: the file exists and is read-only')):
        if os.access(movie, os.W_OK):           # as root
            continue
        capsys.readouterr()
        with pytest.raises(SystemExit) as e:
            main(['nma', pdb, *options, '-o', f'{tmp_path}/new/taken'])
        assert said in str(e.value) and 'analysed' not in capsys.readouterr().out
    locked = tmp_path / 'locked'
    locked.mkdir()
    locked.chmod(0o555)
    try:
        if not os.access(locked, os.W_OK):      # not so as root, or on Windows
            capsys.readouterr()
            with pytest.raises(SystemExit) as e:
                main(['nma', pdb, '-o', f'{locked}/x'])
            assert 'cannot write to' in str(e.value) and 'analysed' not in capsys.readouterr().out
    finally:
        locked.chmod(0o755)


def test_version():
    assert ppitango.__version__


def test_blender_scene(result, tmp_path):
    """The scene opens again with every frame in it, and the contacts drawn change as the frames go by."""
    bpy = pytest.importorskip('bpy')
    pytest.importorskip('molecularnodes')
    from ppitango import blender
    scene = result.scene()
    path = blender.write(scene, str(tmp_path / 'scene.blend'), seconds=2)
    assert [p.name for p in tmp_path.iterdir()] == ['scene.blend']
    import zstandard
    raw = zstandard.ZstdDecompressor().stream_reader(open(path, 'rb')).read()
    assert b'site-packages' not in raw and os.path.expanduser('~').encode() not in raw   # no local paths
    bpy.ops.wm.open_mainfile(filepath=path)
    world = bpy.context.scene
    assert (world.frame_start, world.frame_end) == (1, 2 * blender.FPS) and world.camera is not None
    model = bpy.data.objects['complex']
    assert len(model.data.vertices) == len(scene.ensemble.atoms)
    made = (scene.occ[:, 0] >= 0.5).sum()
    tree = model.modifiers[0].node_group
    drawn = tree.nodes['Contacts']
    for other in [n for n in tree.nodes if n.bl_idname == 'GeometryNodeGroup' and n != drawn]:
        if other.node_tree.name.startswith('Style'):
            tree.nodes.remove(other)                         # leave the contacts alone on the stage
    drawn.inputs['Gap (Å)'].default_value = 1000             # one dash per contact
    world.frame_set(1)
    mesh = model.evaluated_get(bpy.context.evaluated_depsgraph_get()).to_mesh()
    assert len(mesh.vertices) == made * 24                   # a 12-sided cylinder each
