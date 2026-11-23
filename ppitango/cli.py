"""Command line: ppitango nma | md | fep."""
import argparse
import math
import os
import sys
import tempfile

from . import __version__, analysis, ensembles


def _number(kind, least=None, above=None):
    """An argparse type: a number of at least `least`, or greater than `above`."""
    def read(text):
        try:
            value = kind(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f'{text!r} is not a {"whole " if kind is int else ""}number')
        if not math.isfinite(value):        # float() reads "nan" and "inf"
            raise argparse.ArgumentTypeError(f'{text!r} is not a number to work with')
        if (least is not None and value < least) or (above is not None and value <= above):
            raise argparse.ArgumentTypeError(f'must be {f"{least} or more" if least is not None else f"more than {above}"}'
                                             f', got {text}')
        return value
    return read


class _Parser(argparse.ArgumentParser):
    """Says what is wrong with an option in one line, without the list of all options in front."""
    def error(self, message):
        self.exit(2, f'ppitango: {message} (all options: {self.prog} -h)\n')


def _prefix(output, name):
    """Where the outputs go: `output` is a prefix, or a folder (it exists, or ends in a slash) to put files
    called after the input in. The folder is made and tried before any work is done."""
    if output and (output.endswith(('/', os.sep)) or os.path.isdir(output)):
        output = os.path.join(output, name)
    prefix = output or name
    folder = os.path.dirname(os.path.abspath(prefix))
    try:
        os.makedirs(folder, exist_ok=True)
        with tempfile.TemporaryFile(dir=folder):
            pass
    except OSError as e:
        sys.exit(f'ppitango: cannot write to {folder}: {e.strerror}')
    return prefix


def _common(p):
    p.add_argument('--select', nargs=2, metavar=('SIDE_A', 'SIDE_B'),
                   help='the two sides as MDAnalysis selections instead of chain IDs, e.g. '
                        '--select "resid 1-110" "resid 111-199"; a selection without protein in it makes a ligand, '
                        'ion or sugar a side: --select protein "resname LIG"')
    p.add_argument('--renumber', nargs='+', metavar='CHAIN:FIRST', default=[],
                   help='number a chain from this residue number on, e.g. --renumber B:1 (a tpr without its gro, '
                        'and a plain prmtop, count residues straight through all chains)')
    p.add_argument('--chains', nargs=2, metavar=('A', 'B'),
                   help='the two chains whose interface to follow (default: the pair with most contacts; protein against '
                        'nucleic acid when both are there)')
    p.add_argument('-o', '--output', default=None,
                   help='output prefix, or a folder (default: name of the input file, then _nma, _md or _fep; '
                        'for md the input file is the trajectory)')
    out = p.add_argument_group('outputs (default: 2D movie + interactive viewer + CSV tables)')
    out.add_argument('--movie3d', action='store_true',
                     help='also render the split-screen 3D movie: ray-traced models next to the live 2D map (needs PyMOL)')
    out.add_argument('--blender', action='store_true',
                     help='also write a Blender scene (.blend) to restyle and render in Blender: the moving complex, '
                          'its contacts and a camera path (needs Blender as a Python module and Molecular Nodes)')
    out.add_argument('--gif', action='store_true', help='also convert the movies to GIF')
    out.add_argument('--no-movie', action='store_true', help='skip the 2D movie')
    out.add_argument('--no-viewer', action='store_true', help='skip the interactive HTML viewer')
    look = p.add_argument_group('animation')
    look.add_argument('--seconds', type=_number(float, above=0), default=12, help='length of the movies (default: 12)')
    look.add_argument('--fps', type=_number(int, 1), default=30, help='frames per second of the 2D movie (default: 30)')
    look.add_argument('--smooth', type=_number(int, 1), default=None, metavar='FRAMES',
                      help='rolling window that turns flickering contacts into smooth fades '
                           '(default: 5 for MD and FEP, 1 = off for normal modes)')
    look.add_argument('--max-residues', type=_number(int, 1), default=14, help='residues per chain on the 2D map (default: 14)')
    look.add_argument('--title', default='', help='title shown on the movies and the viewer')
    run = p.add_argument_group('analysis')
    run.add_argument('--distance-cutoff', type=_number(float, above=0), default=4.5,
                     help='PandaPPI distance cutoff in Å (default: 4.5)')
    run.add_argument('--workers', type=_number(int, 1), default=None, help='parallel processes (default: all cores)')


def main(argv=None):
    parser = _Parser(
        prog='ppitango', description='ppi-tango: animated protein-protein interface contacts '
                                     'from normal modes, molecular dynamics and free energy perturbation.')
    parser.add_argument('--version', action='version', version=f'ppi-tango {__version__}')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('nma', help='oscillate a complex along an elastic-network normal mode',
                       description='Normal mode analysis needs nothing but a PDB file: the mode is computed here.')
    p.add_argument('pdb', help='PDB file of the complex (or any structure file MDAnalysis reads)')
    p.add_argument('--mode', type=_number(int, 1), default=1, help='which mode, 1 = slowest (default: 1)')
    p.add_argument('--rmsd', type=_number(float, above=0), default=1.0, help='Cα RMSD at full swing, Å (default: 1.0)')
    p.add_argument('--frames', type=_number(int, 2), default=48, help='frames in one oscillation (default: 48)')
    _common(p)

    p = sub.add_parser('md', help='follow the interface along an MD trajectory',
                       description='Any topology and trajectory MDAnalysis reads (pdb, gro, tpr, top, prmtop, psf + xtc, trr, dcd, '
                                   'nc ...), or a multi-model PDB on its own.')
    p.add_argument('topology', help='structure or topology file of the system (pdb, gro, tpr, top, prmtop, psf), '
                                    'or a multi-model PDB on its own. A trajectory saved for part of it (protein '
                                    'only, no water) is matched up')
    p.add_argument('trajectory', nargs='?', help='trajectory file')
    p.add_argument('--waters', action='store_true',
                   help='also follow water-mediated contacts (needs the water in topology and trajectory)')
    p.add_argument('--stride', type=_number(int, 1), default=1, help='read every n-th frame (default: 1)')
    p.add_argument('--max-frames', type=_number(int, 2), default=200, help='analyse at most this many frames (default: 200)')
    p.add_argument('--dt', type=_number(float, above=0), default=None,
                   help='time between saved frames (default: the times stored in the trajectory, else frame numbers)')
    p.add_argument('--unit', default='ns', choices=('fs', 'ps', 'ns', 'us', 'µs', 'ms'),
                   help='unit of --dt: fs, ps, ns, us (or µs), ms (default: ns)')
    _common(p)

    p = sub.add_parser('fep', help='walk through the lambda windows of an alchemical calculation',
                       description='One trajectory (or multi-model PDB) per lambda window, in order.')
    p.add_argument('topology', help='structure or topology file holding the same atoms as the window trajectories')
    p.add_argument('windows', nargs='+', help='one trajectory per lambda window, in order')
    p.add_argument('--waters', action='store_true',
                   help='also follow water-mediated contacts (needs the water in topology and trajectory)')
    p.add_argument('--dg', help='CSV with one row per window and a dG column (cumulative, kcal/mol); '
                                'an optional lambda column sets the axis')
    p.add_argument('--frames-per-window', type=_number(int, 1), default=8, help='frames analysed per window (default: 8)')
    p.add_argument('--label', default='', help='what is being perturbed, e.g. "barstar Y29 → Ala"')
    _common(p)

    args = parser.parse_args(argv)
    source = args.pdb if args.command == 'nma' else args.topology if args.command == 'fep' else \
        args.trajectory or args.topology
    prefix = _prefix(args.output, f'{os.path.splitext(os.path.basename(source))[0]}_{args.command}')
    try:
        try:
            renumber = {c: int(n) for c, n in (item.split(':') for item in args.renumber)} or None
        except ValueError:
            raise ValueError('--renumber takes CHAIN:FIRST items, e.g. --renumber A:3 B:1')
        if args.command == 'nma':
            ens = ensembles.nma(args.pdb, args.chains, args.mode, args.rmsd, args.frames, select=args.select,
                                renumber=renumber)
        elif args.command == 'md':
            ens = ensembles.md(args.topology, args.trajectory, args.chains, args.stride, args.max_frames,
                               args.dt, args.unit, args.select, renumber, args.waters)
        else:
            ens = ensembles.fep(args.topology, args.windows, args.chains, args.dg, args.frames_per_window,
                                args.label, args.select, renumber, args.waters)
    except (ValueError, ImportError, OSError) as e:
        sys.unraisablehook = lambda *_: None   # a reader that could not open its file complains again when dropped
        sys.exit(f'ppitango: {e}')

    try:
        _work(args, ens, prefix)
    except OSError as e:
        sys.exit(f'ppitango: {e.strerror or e}' + (f': {e.filename}' if e.filename else ''))


def _work(args, ens, prefix):
    a, b = ens.chains
    for line in ens.log:
        print(line)
    if ens.n_frames < 2:
        sys.exit('ppitango: only one frame was read, so there is nothing to animate. For a single structure use '
                 'ppitango nma, or PandaPPI itself for a static map.')
    print(f'{ens.title}: {ens.pair}, {ens.n_frames} frames, {len(ens.atoms)} heavy atoms')
    result = analysis.analyze(ens, cutoff=args.distance_cutoff, workers=args.workers)
    if result.contacts.empty:
        sys.exit(f'ppitango: no contacts between {ens.pair} in any frame - nothing to animate.')
    scene = result.scene(smooth=args.smooth, max_residues=args.max_residues)
    files = result.save(prefix)

    if not args.no_viewer:
        from . import viewer
        files.append(viewer.write(scene, f'{prefix}_viewer.html', args.title))
    from . import movie
    movies = []
    if not args.no_movie:
        print('  rendering 2D movie ...')
        movies.append(movie.render(scene, f'{prefix}_contacts.mp4', args.seconds, args.fps, args.title))
    if args.movie3d:
        from . import movie3d
        try:
            movies.append(movie3d.render(scene, f'{prefix}_3d.mp4', args.seconds, title=args.title))
        except RuntimeError as e:
            print(f'ppitango: 3D movie skipped: {e}', file=sys.stderr)
    files += movies + ([movie.to_gif(m) for m in movies] if args.gif else [])
    if args.blender:
        from . import blender
        try:
            files.append(blender.write(scene, f'{prefix}.blend', args.seconds))
        except RuntimeError as e:
            print(f'ppitango: Blender scene skipped: {e}', file=sys.stderr)
    print('Wrote:')
    for f in files:
        print(f'  {f}')


if __name__ == '__main__':
    main()
