# ppi-tango

**It takes two.** Animated protein-protein interface contacts from normal modes, molecular
dynamics and free energy perturbation.

[PandaPPI](https://github.com/pritampanda15/PandaPPI) draws every contact at an interface from one
static structure. ppi-tango runs it on every frame of a moving complex and turns the result into
things you can put in a talk:

![Split-screen 3D movie of barnase-barstar](demos/barnase_barstar_md_3d.gif)

| Output | File | What it is |
| --- | --- | --- |
| **2D movie** | `*_contacts.mp4` (`.gif`) | PandaPPI's interaction map, alive: contact lines fade in and out, distances tick, bars and timelines follow along |
| **Interactive 3D viewer** | `*_viewer.html` | The real models moving, every contact drawn between its two atoms with a live distance, camera presets, a guided tour, and the 2D map in sync. One self-contained file, works offline |
| **3D movie** | `*_3d.mp4` | Ray-traced models next to the live 2D map; the camera flies in to the interface and around it. Needs PyMOL |
| **Tables** | `*_contacts.csv`, `*_buried_surface.csv` | Every contact and every buried residue in every frame |

## Install

```bash
pip install "ppi-tango[md] @ git+https://github.com/varunsharm16/ppi-tango"
```

That is everything for normal modes, MD and FEP input, the 2D movie and the interactive viewer
(ffmpeg comes bundled). Leave out `[md]` if you only need normal modes or multi-model PDB files.
The 3D movie additionally needs the `pymol` command:

```bash
conda install -c conda-forge pymol-open-source
```

## Use

### Normal mode analysis: just a PDB file

```bash
ppitango nma 1BRS.pdb --chains A D
```

The elastic network model is built in: no other program needed. `--mode 2` picks another mode,
`--rmsd 1.5` a larger swing.

### Molecular dynamics: a topology and a trajectory

```bash
ppitango md top.pdb traj.dcd --chains A D --dt 0.1 --unit ns
```

`top.pdb` must carry chain IDs and match the trajectory's atoms (water, ions and hydrogens are
dropped automatically). Anything [MDAnalysis](https://www.mdanalysis.org) reads works: xtc, dcd,
trr, nc. A multi-model PDB on its own works too. A complex split across the periodic box is put
back together, and frames are fitted onto the first.

### Free energy perturbation: one trajectory per λ window

```bash
ppitango fep top.pdb window_*.dcd --chains A D --dg dg.csv --label "barstar Y29 side chain off"
```

`dg.csv` has one row per window and a `dG` column (cumulative free energy, kcal/mol, from
alchemlyb / pymbar / your FEP package). The free energy curve is drawn against λ next to the
contacts that go with it.

### Options that matter

| Option | |
| --- | --- |
| `--chains A B` | the two chains to follow (default: the pair with most contacts) |
| `--movie3d` | also render the split-screen 3D movie |
| `--gif` | also write the 2D movie as a GIF |
| `--smooth N` | rolling window over N frames. Contacts hovering around a cutoff would flicker; this turns flicker into fades (default 5 for MD/FEP, off for normal modes) |
| `--seconds`, `--fps` | movie length and frame rate |
| `--max-residues` | residues per chain on the 2D map (default 14, the busiest ones) |
| `--title` | title on the movies and viewer |

`ppitango nma --help` (and `md`, `fep`) lists the rest.

### In slides

The movies are H.264 MP4 and play in PowerPoint, Keynote and Google Slides. The viewer is a
single HTML file, so in any web-based deck (reveal.js, including
[live slides](https://github.com/ramithuh/protein-design-interactive-slides)) it is one line:

```html
<iframe src="barnase_barstar_md_viewer.html#autoplay" style="width:100%;height:100%;border:0"></iframe>
```

`#autoplay` starts playback, `#tour` also starts the camera tour. In the viewer: drag to rotate,
scroll to zoom, click a contact or residue on the map to fly to it, click a legend chip to hide a
contact type, space to play, arrow keys to step.

### From Python

```python
import ppitango

ens = ppitango.md('top.pdb', 'traj.dcd', chains=['A', 'D'], dt=0.1)   # or ppitango.nma(...), ppitango.fep(...)
result = ppitango.analyze(ens)            # PandaPPI on every frame, in parallel
result.contacts                           # tidy DataFrame: frame, x, type, residue_a, atom_a, residue_b, atom_b, distance
scene = result.scene(smooth=5)            # per-contact occupancy and distance over frames

from ppitango import movie, viewer, movie3d
movie.render(scene, 'contacts.mp4')
viewer.write(scene, 'viewer.html')
movie3d.render(scene, '3d.mp4')
```

## Demos

Everything in [`demos/`](demos) was made with the commands in [`examples/make_demos.sh`](examples/make_demos.sh).
Open the `*_viewer.html` files in a browser after cloning (GitHub does not render them in place).

| Complex | Motion | Files |
| --- | --- | --- |
| Barnase-barstar (1BRS, chains A/D) | 20 ns explicit-solvent MD | `barnase_barstar_md_*` |
| Barnase-barstar | FEP: barstar Tyr29 side chain switched off over 12 λ windows | `barnase_barstar_fep_*` |
| Barnase-barstar | normal mode 1 | `barnase_barstar_nma_*` |
| Fab-VEGF (1BJ1, chains H/W), PandaPPI's own example | normal mode 1 | `fab_vegf_nma_*` |
| ACE2-SARS-CoV-2 RBD (6M0J, chains A/E) | normal mode 1 | `ace2_rbd_nma_*` |

![2D movie of barnase-barstar along its slowest normal mode](demos/barnase_barstar_nma_contacts.gif)

The MD and FEP runs are in [`examples/simulations/simulate.py`](examples/simulations/simulate.py)
(OpenMM, Amber14SB/TIP3P, 150 mM NaCl, 298 K). They are short demonstration runs, not
converged production data: see the notes below.

## How it works, and what to keep in mind

- **Contacts** are PandaPPI's, with its definitions and cutoffs, computed per frame on the heavy
  atoms of the two chains. Hydrogen bonds are therefore distance-based (donor-acceptor ≤ 3.5 Å),
  and water-mediated contacts are not tracked because frames are stripped of solvent. Van der
  Waals atom contacts are counted but not drawn, as in PandaPPI.
- **Occupancy.** A contact is one (type, residue, residue) triple. Its presence over frames is
  smoothed with a rolling mean into an occupancy between 0 and 1, which sets the opacity of its
  line. Its distance is measured in every frame between the two atoms that make it, so the number
  keeps ticking while the line fades.
- **Normal modes** come from an anisotropic network model on Cα atoms (15 Å cutoff). Each residue
  moves rigidly with its Cα along the mode, which is a linear extrapolation: fine for showing which
  contacts a collective motion strains, wrong for large amplitudes (bonds between neighbouring
  residues stretch, atoms can overlap). Keep `--rmsd` near 1 Å. Amplitudes and time scales of a
  mode are not physical.
- **FEP demo.** The example calculation annihilates the Tyr29 side chain beyond Cβ (a
  pseudo-alanine) with 1 ns per window in complex and free barstar, analysed with MBAR. It is
  there to show the animation, and its ΔΔG is a rough estimate: sampling is short, and a proper
  mutation would use a hybrid topology.
- **One chain pair at a time.** For an antibody, run heavy-antigen and light-antigen separately.
- **Residue numbers** follow PandaPPI (name + number, insertion codes ignored).

## Credits

Contact detection and the visual language (colours, badges, the map itself) are
[PandaPPI](https://github.com/pritampanda15/PandaPPI)'s, by Pritam Kumar Panda. The viewer uses
[3Dmol.js](https://3dmol.csb.pitt.edu) (BSD-3-Clause, bundled in `ppitango/vendor`); the 3D movie
uses [PyMOL](https://github.com/schrodinger/pymol-open-source). MIT licence.
