<p align="center"><img src="docs/header.webp" alt="ppi-tango: it takes two" width="100%"></p>

# ppi-tango

**It takes two.** Animated protein-protein interface contacts from normal modes, molecular
dynamics and free energy perturbation.

[PandaPPI](https://github.com/pritampanda15/PandaPPI) draws every contact at an interface from one
static structure. ppi-tango runs it on every frame of a moving complex and turns the result into
things you can put in a talk:

![Split-screen 3D movie of barnase-barstar over 20 ns of molecular dynamics](docs/hero.gif)

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
The 3D movie additionally needs PyMOL, from pip or from conda:

```bash
pip install "ppi-tango[md,movie3d] @ git+https://github.com/varunsharm16/ppi-tango"
conda install -c conda-forge pymol-open-source      # or this, next to a plain install
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
ppitango md top.pdb traj.dcd --chains A D --dt 0.1          # OpenMM, NAMD: PDB + DCD
ppitango md md.tpr md.xtc                                   # GROMACS (or md.gro, or topol.top, + md.xtc)
ppitango md complex.prmtop md.nc                            # Amber
ppitango md step5.psf step5.dcd --dt 0.1                    # CHARMM / NAMD
```

Give the files as your MD program wrote them. The first file is any structure or topology
[MDAnalysis](https://www.mdanalysis.org) reads (pdb, gro, tpr, top, prmtop/parm7, psf); the second any
trajectory (xtc, trr, dcd, nc ...). Nothing needs preparing by hand:

- **Water, ions, lipids, ligands and hydrogens** are dropped. What stays is every residue with a protein or
  nucleic acid backbone, whatever the force field calls it (HIE, HSD, CYX ... are understood).
- **Periodic boxes**: chains cut by the box are made whole and the partners are put back next to each
  other in every frame, so a raw trajectory works without `gmx trjconv -pbc` or `cpptraj autoimage`.
- **Topology and trajectory for different sets of atoms.** A trajectory saved for part of the system
  (`gmx trjconv` or `compressed-x-grps` with Protein, `cpptraj strip`) works with the full topology: the part
  is recognised by its atom count (protein, protein and nucleic acid, everything but water, everything but
  water and ions, each also without hydrogens, the backbone, the Cα atoms) and checked against the bond
  lengths. The reverse works too: a solute-only topology with a full trajectory, when the solute comes first
  in the system, as MD programs write it. What was assumed is printed.
- **Chains.** With chain IDs in the topology, `--chains A D` picks two; left out, the pair with most contacts
  is followed. Files without chain IDs (gro, prmtop, most psf) get their chains lettered A, B, C ... in file
  order, split wherever the backbone breaks, and the table is printed so you can pick:

  ```
  Chains (no chain IDs in the topology: chains lettered in file order, split where the backbone breaks):
    A: 108 residues, VAL3 to ARG110
    B: 87 residues, LYS1 to SER89
  ```

- **Anything else is a selection**: `--select "segid PROA" "segid PROB"`, `--select "resid 1-220" "resid 221-330"`,
  or a side made of several chains, `--select "chainID H L" "chainID W"` for an antibody against its antigen
  ([selection language](https://userguide.mdanalysis.org/stable/selections.html)).
- **Residue numbers** are the ones of your structure wherever a file still has them: PDB, `.gro`, a GROMACS
  `.top` (each molecule keeps its numbers), a `.tpr` with the `.gro` of the same name next to it, an Amber
  `prmtop` that carries them (`parmed`'s `addPDB`). A `.tpr` or `prmtop` on its own counts 1, 2, 3 ...
  straight through all chains: `--renumber B:1` (chain, number of its first residue) puts a chain back on its
  own numbers. Insertion codes are kept (Tyr 100A), and when a side holds several chains its residues are
  named with the chain in front (H:100A, L:32).
- **Ligands, ions, sugars, lipids**: name them in a selection of their own and they are a side like any
  other: `--select protein "resname JZ4"`, `--select "chainID A" "resname ZN"`, `--select protein "resname NAG"`.
  A selection that holds protein or nucleic acid keeps only that (so `chainID A` does not drag the waters
  and buffer of chain A along).
- **Water-mediated contacts**: `--waters` follows the waters within hydrogen-bond reach of both sides in every
  frame (give the full-system topology and trajectory, not a stripped one).
- **Protein-DNA and protein-RNA**: a file with both is read as all protein chains against all nucleic acid
  chains unless you say otherwise. Hydrogen bonds come from PandaPPI; ppi-tango adds salt bridges from Arg,
  Lys and His to backbone phosphates, stacking of aromatic side chains on bases and cation-π from Arg and Lys
  to bases. Works the same for `nma` (the network gets a node at C4' of every nucleotide) and `fep`.
- **Time axis**: read from the trajectory when it stores times (xtc, trr, nc). DCD files do not, so give
  `--dt` (time between saved frames, with `--unit`), or the axis counts frames.
- **Long trajectories**: 200 evenly spaced frames are analysed (`--max-frames`, `--stride`), and only those
  are read, so trajectory length hardly matters.

A multi-model PDB on its own works too.

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
| `--gif` | also convert the movies to GIF (for READMEs and chat; use the MP4 in slides) |
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

### As a talk you walk through with the arrow keys

[`slides/talks/barnase_barstar`](slides/talks/barnase_barstar) is a finished ten-slide talk built from the MD
and FEP demos: open `slides/talks/barnase_barstar/index.html` in a browser and press → until it ends. Every
press plays one step: [manim-web](https://github.com/maloyan/manim-web) animations for the diagrams, the real
cartoons following the trajectory in 3D with distances measured live, a camera that moves into the interface.
Nothing else is needed while presenting; the 3D slides can still be dragged and zoomed if someone asks.

![One aspartate, four partners](docs/talk_hub.png)

It is built on Ramith Hettiarachchi's
[live slides template](https://github.com/ramithuh/protein-design-interactive-slides) (a copy is in `slides/`;
[`slides/AGENTS.md`](slides/AGENTS.md) has its conventions). To make the same kind of talk from your own run:

```bash
python slides/tools/prep_interface.py top.pdb traj.dcd --chains A D --dt 0.1 --focus D:39 D:29 \
    --dg dg.csv slides/talks/my_talk/data/interface.json
```

writes the one data file the slides draw (cartoons, the trajectory, every contact with its occupancy, the
distances of the `--focus` residues' partners, the FEP curve; needs `mdtraj`). It takes the same input files
as `ppitango md`. Copy
`slides/talks/barnase_barstar` to `slides/talks/my_talk`, keep `slides/stage3d.js` (the 3D stage) and rewrite
the story in the other slide files: the numbers on them are read from the data, the sentences are yours.
`python slides/tools/check_deck.py slides/talks/my_talk` steps through every slide in a headless browser.

To send a talk to someone, pack it into one file that needs nothing next to it and no network (about 7 MB):

```bash
python slides/tools/bundle.py slides/talks/barnase_barstar
```

writes `slides/talks/barnase_barstar/barnase_barstar.html`: attach it to a message or put it on any web host.
Whoever opens it clicks the slide once, then uses the arrow keys.

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
(OpenMM, Amber14SB/TIP3P, 150 mM NaCl, 298 K); their trajectories are in
`examples/simulations/output`, so `examples/make_demos.sh` runs as is. They are short
demonstration runs, not converged production data: see the notes below.

## What was tested, and known limits

Version 0.1 was run on the same 195-residue complex as four MD programs write it, and each input had to give
the contacts of the reference run (`tests/test_trajectory.py` keeps the core of this, and runs on every push):

| Input | Result |
| --- | --- |
| GROMACS 2025: `md.tpr` or `md.gro` + raw `md.xtc` (dodecahedron box, water, ions, molecules cut by the box) | same contacts and distances as the `trjconv -pbc mol` trajectory |
| GROMACS 2026, with and without GROMACS installed, MDAnalysis 2.10 and 2.11: `md.tpr`, `md.gro` or `topol.top` (with the `.itp` files `pdb2gmx` wrote) + raw `md.xtc` | same contacts, distances and residue numbers from all three |
| The same run with the trajectory saved for Protein, for Protein-H or without water, against the full `.tpr`, `.gro` and `.top`; and a protein-only `.gro` with the full trajectory | same again, all ten combinations |
| GROMACS protein-only `.gro`/`.pdb` + `.xtc` | same |
| Amber: `prmtop` + NetCDF, full system (truncated octahedron) and stripped; `ambpdb`-style PDB without chains | all four agree |
| CHARMM-style PSF + DCD; PDB with HSD, ILE CD, OT1/OT2 and segment names only | same as reference |
| OpenMM PDB + DCD, multi-model PDB, FEP windows | unchanged from the demos |
| 12,800-frame trajectory (150 MB) | 70 s on 8 cores without movies |
| Antibody heavy + light chain against antigen (`--select`, 1BJ1); Fab with Kabat insertion codes (1HZH, 100A-100J) | run; residues keep their names (H:100A) |
| Ligand in a buried pocket (T4 lysozyme L99A with 2-propylphenol, 3HTB); zinc sites of a zinc finger (1ZAA) and of ACE2 (6M0J) | hydrogen bonds and 13 hydrophobic contacts to the ligand; every zinc on its four Cys/His |
| Protein-DNA: zinc finger Zif268 on its site (1ZAA), λ repressor-operator (1LMB) with 50 ps of explicit-solvent MD, normal modes, and the FEP path on stand-in windows | run; phosphate salt bridges, base stacking, cation-π and water-mediated contacts found |
| Python 3.9 and 3.12 on macOS; 3.10, 3.12, 3.13 on Linux, 3.12 on macOS and Windows (CI) | pass |
| Alliance Canada cluster (narval) login node, no display: fresh virtualenv from the site wheels (numpy 2.4, MDAnalysis 2.10); GROMACS 2025 and 2026 `.tpr`, `.gro`, `topol.top` with and without the GROMACS module loaded, full and protein-only trajectories, `--waters`, normal modes, 2D movie | run; all nine topology/trajectory combinations give the same contacts and residue numbers |

What does not work, or needs care:

- **A trajectory for an unusual part of the system.** A trajectory and a topology with different atom
  counts are matched up as described above. A custom index group (one chain of several, a hand-picked
  region) is not recognised and stops with an error saying both counts: write the topology for that group
  (`gmx editconf -n`, `gmx convert-tpr -n`).
- **A GROMACS `.top` needs a trajectory**, since it holds no coordinates (normal modes need the `.gro`). Its
  `#include` lines are followed. Without GROMACS on the machine the force field files are not there: water and
  ions are then filled in from their names, which covers the usual protein-in-water system; another molecule
  defined only in a force field file stops with an error naming it (load GROMACS, or give the `.gro`).
- **A `.tpr` from a very new GROMACS** may be ahead of MDAnalysis: GROMACS 2026 files are read from
  MDAnalysis 2.11 on, not by 2.10; GROMACS 2025 files not by 2.7. ppi-tango then uses the `.gro` or `.top`
  with the same name next to the `.tpr` (or `topol.top`) if there is one, and says so; otherwise write a
  `.gro` (`gmx editconf -f md.tpr -o md.gro`) and give that, or `pip install -U MDAnalysis`.
- **Residue numbers of a lone `.tpr` or a plain `prmtop`** run straight through the chains (see above):
  keep the `.gro` next to the `.tpr`, or use `--renumber`.
- **xtc rounds coordinates to 0.01 Å.** A contact sitting exactly on a cutoff can differ between an xtc and
  a full-precision file of the same run; occupancies are not affected in any way that shows.
- **Ligands, ions, sugars and lipids are a side only on their own.** A side cannot mix protein with its
  ligand or its glycans. For such a side ppi-tango reports PandaPPI's hydrogen bonds (by distance between
  N and O atoms, without hydrogens), plus its own hydrophobic contacts (carbons without N or O attached, to
  side-chain carbons of hydrophobic residues), metal coordination and halide salt bridges. Not detected for
  ligands: aromatic stacking, cation-π, halogen bonds, covalent links. Their atom names are taken as they
  are, so two ligand atoms with one name in one residue cannot be told apart.
- **Water-mediated contacts** need `--waters` and the water in both files (not possible when a solute-only
  topology is matched to a full trajectory). The bridging water is drawn in the 3D viewer and movie; the
  distance shown is still the one between the two residues' atoms.
- **Protein-nucleic acid** complexes were tested on structures and a short MD run only. No alchemical
  protein-DNA calculation was run: the FEP path was checked for reading and drawing with stand-in windows.
- **One frame is not an animation**: `md` on a single structure stops; use `nma`, or PandaPPI for a static map.
- **More than 9999 residues on one side** cannot be written to the PDB frames PandaPPI reads.
- **3D movie** needs PyMOL and is the slow part (about 3 minutes for a 12 s movie on 8 cores); everything
  else runs headless on a cluster login node without it.
- **Windows**: the tests, the README commands and a 3D movie (PyMOL from pip) run in CI on every push; nobody
  has looked at the pictures made there.
- **Slides data** (`slides/tools/prep_interface.py`) reads the same inputs (`--select`, `--renumber` too).
  The deck was made for two proteins: a nucleic acid side is written as a plain backbone trace through C4',
  and the data file was produced for a protein-DNA run but the deck was not viewed with it. A ligand or ion
  side is not supported there.

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
- **FEP demo.** The example calculation annihilates the barstar Tyr29 side chain beyond Cβ (a
  pseudo-alanine) with 1 ns per window in the complex and in free barstar, analysed with MBAR. It
  gives ΔΔG = 4.0 ± 0.6 kcal/mol; the measured value for Y29A is about 3.4 kcal/mol (Schreiber &
  Fersht, 1995). Treat the agreement as encouraging, not as validation: sampling is short, the
  error bar is statistical only, and a proper mutation would use a hybrid topology.
- **FEP contacts are geometric.** Alchemical atoms stay in the coordinates while they are being
  switched off, so a contact drawn at high λ means "these atoms still sit close", not "this
  interaction still contributes". Watch what the neighbours do as λ grows.
- **One chain pair at a time.** For an antibody, run heavy-antigen and light-antigen separately.
- **Residue numbers** follow PandaPPI (name + number, insertion codes ignored).

## Citing

If ppi-tango helps your work, cite this repository and PandaPPI, which does the contact detection:

```
Varun Sharma. ppi-tango: animated protein-protein interface contacts from normal modes, molecular
dynamics and free energy perturbation. GitHub repository, https://github.com/varunsharm16/ppi-tango

Pritam Kumar Panda. PandaPPI: mapping and visualizing molecular interactions at protein-protein,
protein-nucleic acid and antigen-antibody interfaces. GitHub repository,
https://github.com/pritampanda15/PandaPPI
```

### References

Structures used in the demos (coordinates from the [RCSB PDB](https://www.rcsb.org)):

- **1BRS**, barnase-barstar: Buckle AM, Schreiber G, Fersht AR (1994). Protein-protein recognition: crystal
  structural analysis of a barnase-barstar complex at 2.0-Å resolution. *Biochemistry* 33, 8878-8889.
- **1BJ1**, Fab-VEGF: Muller YA, Chen Y, Christinger HW, Li B, Cunningham BC, Lowman HB, de Vos AM (1998). VEGF
  and the Fab fragment of a humanized neutralizing antibody: crystal structure of the complex at 2.4 Å
  resolution and mutational analysis of the interface. *Structure* 6, 1153-1167.
- **6M0J**, ACE2-SARS-CoV-2 RBD: Lan J, Ge J, Yu J, et al. (2020). Structure of the SARS-CoV-2 spike
  receptor-binding domain bound to the ACE2 receptor. *Nature* 581, 215-220.

Barnase-barstar numbers quoted in the talk and in the notes above:

- Affinity: Schreiber G, Fersht AR (1993). Interaction of barnase with its polypeptide inhibitor barstar
  studied by protein engineering. *Biochemistry* 32, 5145-5150. Hartley RW (1993). Directed mutagenesis and
  barnase-barstar recognition. *Biochemistry* 32, 5978-5984.
- Y29A and the other interface mutants: Schreiber G, Fersht AR (1995). Energetics of protein-protein
  interactions: analysis of the barnase-barstar interface by single mutations and double mutant cycles.
  *J Mol Biol* 248, 478-486.

Methods:

- Normal modes (anisotropic network model): Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, Bahar I
  (2001). Anisotropy of fluctuation dynamics of proteins with an elastic network model. *Biophys J* 80, 505-515.
- Simulations: Eastman P, et al. (2024). OpenMM 8: molecular dynamics simulation with machine learning
  potentials. *J Phys Chem B* 128, 109-116.
- Force field and water: Maier JA, et al. (2015). ff14SB: improving the accuracy of protein side chain and
  backbone parameters from ff99SB. *J Chem Theory Comput* 11, 3696-3713. Jorgensen WL, et al. (1983). Comparison
  of simple potential functions for simulating liquid water. *J Chem Phys* 79, 926-935.
- Free energies (MBAR): Shirts MR, Chodera JD (2008). Statistically optimal analysis of samples from multiple
  equilibrium states. *J Chem Phys* 129, 124105.

Software ppi-tango reads with or draws with:

- MDAnalysis: Michaud-Agrawal N, Denning EJ, Woolf TB, Beckstein O (2011). MDAnalysis: a toolkit for the
  analysis of molecular dynamics simulations. *J Comput Chem* 32, 2319-2327.
- MDTraj (slide data, secondary structure): McGibbon RT, et al. (2015). MDTraj: a modern open library for the
  analysis of molecular dynamics trajectories. *Biophys J* 109, 1528-1532. Kabsch W, Sander C (1983).
  Dictionary of protein secondary structure. *Biopolymers* 22, 2577-2637.
- 3Dmol.js: Rego N, Koes D (2015). 3Dmol.js: molecular visualization with WebGL. *Bioinformatics* 31, 1322-1324.
- PyMOL: The PyMOL Molecular Graphics System, Schrödinger, LLC.
- Matplotlib: Hunter JD (2007). Matplotlib: a 2D graphics environment. *Comput Sci Eng* 9, 90-95.

## Credits

Contact detection and the visual language (colours, badges, the map itself) are
[PandaPPI](https://github.com/pritampanda15/PandaPPI)'s, by Pritam Kumar Panda. The viewer uses
[3Dmol.js](https://3dmol.csb.pitt.edu) (BSD-3-Clause, bundled in `ppitango/vendor`); the 3D movie
uses [PyMOL](https://github.com/schrodinger/pymol-open-source). The talk in `slides/` is built on the
[live slides template](https://github.com/ramithuh/protein-design-interactive-slides) (MIT), which bundles manim-web,
three.js, Reveal.js and the Lato font with their licences in `slides/vendor`. MIT licence.
