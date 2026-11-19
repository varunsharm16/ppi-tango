"""Command line: ppitango nma | md | fep."""
import argparse
import os
import sys

from . import __version__, analysis, ensembles


def _common(p):
    p.add_argument('--select', nargs=2, metavar=('SIDE_A', 'SIDE_B'),
                   help='the two sides as MDAnalysis selections instead of chain IDs, e.g. '
                        '--select "resid 1-110" "resid 111-199"')
    p.add_argument('--renumber', nargs='+', metavar='CHAIN:FIRST', default=[],
                   help='number a chain from this residue number on, e.g. --renumber B:1 (tpr and prmtop files count '
                        'residues straight through all chains)')
    p.add_argument('--chains', nargs=2, metavar=('A', 'B'),
                   help='the two chains whose interface to follow (default: the pair with most contacts)')
    p.add_argument('-o', '--output', default=None, help='output prefix (default: name of the input file)')
    out = p.add_argument_group('outputs (default: 2D movie + interactive viewer + CSV tables)')
    out.add_argument('--movie3d', action='store_true',
                     help='also render the split-screen 3D movie: ray-traced models next to the live 2D map (needs PyMOL)')
    out.add_argument('--gif', action='store_true', help='also convert the movies to GIF')
    out.add_argument('--no-movie', action='store_true', help='skip the 2D movie')
    out.add_argument('--no-viewer', action='store_true', help='skip the interactive HTML viewer')
    look = p.add_argument_group('animation')
    look.add_argument('--seconds', type=float, default=12, help='length of the movies (default: 12)')
    look.add_argument('--fps', type=int, default=30, help='frames per second of the 2D movie (default: 30)')
    look.add_argument('--smooth', type=int, default=None, metavar='FRAMES',
                      help='rolling window that turns flickering contacts into smooth fades '
                           '(default: 5 for MD and FEP, 1 = off for normal modes)')
    look.add_argument('--max-residues', type=int, default=14, help='residues per chain on the 2D map (default: 14)')
    look.add_argument('--title', default='', help='title shown on the movies and the viewer')
    run = p.add_argument_group('analysis')
    run.add_argument('--distance-cutoff', type=float, default=4.5,
                     help='PandaPPI distance cutoff in Å (default: 4.5)')
    run.add_argument('--workers', type=int, default=None, help='parallel processes (default: all cores)')


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='ppitango', description='ppi-tango: animated protein-protein interface contacts '
                                     'from normal modes, molecular dynamics and free energy perturbation.')
    parser.add_argument('--version', action='version', version=f'ppi-tango {__version__}')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('nma', help='oscillate a complex along an elastic-network normal mode',
                       description='Normal mode analysis needs nothing but a PDB file: the mode is computed here.')
    p.add_argument('pdb', help='PDB file of the complex (or any structure file MDAnalysis reads)')
    p.add_argument('--mode', type=int, default=1, help='which mode, 1 = slowest (default: 1)')
    p.add_argument('--rmsd', type=float, default=1.0, help='Cα RMSD at full swing, Å (default: 1.0)')
    p.add_argument('--frames', type=int, default=48, help='frames in one oscillation (default: 48)')
    _common(p)

    p = sub.add_parser('md', help='follow the interface along an MD trajectory',
                       description='Any topology and trajectory MDAnalysis reads (pdb, gro, tpr, top, prmtop, psf + xtc, trr, dcd, '
                                   'nc ...), or a multi-model PDB on its own.')
    p.add_argument('topology', help='structure or topology file holding the same atoms as the trajectory '
                                    '(pdb, gro, tpr, top, prmtop, psf), or a multi-model PDB on its own')
    p.add_argument('trajectory', nargs='?', help='trajectory file')
    p.add_argument('--waters', action='store_true',
                   help='also follow water-mediated contacts (needs the water in topology and trajectory)')
    p.add_argument('--stride', type=int, default=1, help='read every n-th frame (default: 1)')
    p.add_argument('--max-frames', type=int, default=200, help='analyse at most this many frames (default: 200)')
    p.add_argument('--dt', type=float, default=None,
                   help='time between saved frames (default: the times stored in the trajectory, else frame numbers)')
    p.add_argument('--unit', default='ns', help='unit of --dt (default: ns)')
    _common(p)

    p = sub.add_parser('fep', help='walk through the lambda windows of an alchemical calculation',
                       description='One trajectory (or multi-model PDB) per lambda window, in order.')
    p.add_argument('topology', help='structure or topology file holding the same atoms as the window trajectories')
    p.add_argument('windows', nargs='+', help='one trajectory per lambda window, in order')
    p.add_argument('--waters', action='store_true',
                   help='also follow water-mediated contacts (needs the water in topology and trajectory)')
    p.add_argument('--dg', help='CSV with one row per window and a dG column (cumulative, kcal/mol); '
                                'an optional lambda column sets the axis')
    p.add_argument('--frames-per-window', type=int, default=8, help='frames analysed per window (default: 8)')
    p.add_argument('--label', default='', help='what is being perturbed, e.g. "barstar Y29 → Ala"')
    _common(p)

    args = parser.parse_args(argv)
    try:
        try:
            renumber = {c: int(n) for c, n in (item.split(':') for item in args.renumber)} or None
        except ValueError:
            raise ValueError('--renumber takes CHAIN:FIRST items, e.g. --renumber A:3 B:1')
        if args.command == 'nma':
            ens = ensembles.nma(args.pdb, args.chains, args.mode, args.rmsd, args.frames, select=args.select,
                                renumber=renumber)
            source = args.pdb
        elif args.command == 'md':
            ens = ensembles.md(args.topology, args.trajectory, args.chains, args.stride, args.max_frames,
                               args.dt, args.unit, args.select, renumber, args.waters)
            source = args.trajectory or args.topology
        else:
            ens = ensembles.fep(args.topology, args.windows, args.chains, args.dg, args.frames_per_window,
                                args.label, args.select, renumber, args.waters)
            source = args.topology
    except (ValueError, ImportError, FileNotFoundError) as e:
        sys.exit(f'ppitango: {e}')

    prefix = args.output or f'{os.path.splitext(os.path.basename(source))[0]}_{args.command}'
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    a, b = ens.chains
    for line in ens.log:
        print(line)
    if ens.n_frames < 2:
        sys.exit('ppitango: only one frame was read, so there is nothing to animate. For a single structure use '
                 'ppitango nma, or PandaPPI itself for a static map.')
    print(f'{ens.title}: chains {a} and {b}, {ens.n_frames} frames, {len(ens.atoms)} heavy atoms')
    result = analysis.analyze(ens, cutoff=args.distance_cutoff, workers=args.workers)
    if result.contacts.empty:
        sys.exit(f'ppitango: no contacts between chains {a} and {b} in any frame - nothing to animate.')
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
    print('Wrote:')
    for f in files:
        print(f'  {f}')


if __name__ == '__main__':
    main()
