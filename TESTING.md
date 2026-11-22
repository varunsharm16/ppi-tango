# Testing ppi-tango from the outside

For anyone (person or agent) who wants to check that ppi-tango does what its [README](README.md) says. You
play a user who has MD data and wants the movies. Everything needed is in this repository: no downloads
besides the install, no account, no cluster.

The goal is to find where the tool or its documentation is wrong. A command that fails, an output that does
not match the README, a sentence that promises more than the tool does: all of these are findings.

## 1. Install, as the README says

Use a fresh environment each time, so that nothing already on the machine helps.

```bash
git clone https://github.com/varunsharm16/ppi-tango && cd ppi-tango
python -m venv .venv && source .venv/bin/activate          # Python 3.10 to 3.13
pip install ".[md,test,movie3d]"
ppitango --version                                          # ppi-tango 0.1.0
pytest -q                                                   # all pass; one Blender test is skipped
```

The README's own install line (`pip install "ppi-tango[md] @ git+https://github.com/varunsharm16/ppi-tango"`)
must work too, without the clone. The Blender scene needs a second environment with Python 3.11 exactly
(section 4).

## 2. The README commands on the shipped runs

Barnase-barstar (PDB 1BRS, chains A and D). The MD and FEP runs are short demonstration runs.

```bash
ppitango nma examples/data/1BRS.pdb --chains A D -o out/nma
ppitango md  examples/simulations/output/md_1brs/top.pdb examples/simulations/output/md_1brs/traj.dcd \
    --chains A D --dt 0.1 --movie3d --gif -o out/md
ppitango fep examples/simulations/output/fep_1brs/top.pdb examples/simulations/output/fep_1brs/window_*.dcd \
    --chains A D --dg examples/simulations/output/fep_1brs/dg.csv --frames-per-window 10 -o out/fep
```

What we get, to compare with:

| Run | Line printed | Rows in `*_contacts.csv` |
| --- | --- | --- |
| nma | `1BRS · normal mode 1: chains A and D, 48 frames, 1557 heavy atoms` | 2425 |
| md | `traj · molecular dynamics: chains A and D, 200 frames, 1565 heavy atoms` | 9813 |
| fep | 12 windows, 10 frames each | 5929 |

Each run writes `*_contacts.mp4`, `*_viewer.html`, `*_contacts.csv` and `*_buried_surface.csv`; `--movie3d`
adds `*_3d.mp4`, `--gif` the GIFs. Our own outputs of these commands are in [`demos/`](demos) to compare by eye.
Things to check in them:

- the movies play, and the contact lines on the 2D map fade in and out rather than flicker;
- the viewer opens offline in a browser, the models move, each contact is drawn between two atoms (to the
  centre of the ring for cation-π, CH-π, π-π and amide-aromatic contacts, whose `atom_a` or `atom_b` is
  empty in the CSV), and the browser console shows no error;
- a contact drawn in a frame is in the CSV for that frame, with that distance;
- residue names and numbers are those of the PDB file (barnase Arg 59, Arg 83, Arg 87, His 102; barstar
  Asp 35, Asp 39, Glu 76 are the well-known interface residues and should be there).

## 3. One system, saved the way each MD program saves it

A user's files are rarely a tidy PDB. [`examples/formats/`](examples/formats) holds one short GROMACS run of
the same complex in water (37044 atoms, dodecahedron box, 11 frames over 20 ps, GROMACS 2026, Amber99SB-ILDN,
TIP3P), and that run again as other programs and habits leave it:

| Files | What makes it hard |
| --- | --- |
| `gromacs/md.gro`, `md.tpr`, `topol.top` (+ two `.itp`) with `md.xtc` | no chain IDs; water and ions; both chains cut by the box in every frame; the `.tpr` is newer than MDAnalysis before 2.11 reads; the `.top` includes force field files you may not have |
| `gromacs/md_whole.xtc` | the same trajectory after `gmx trjconv -pbc mol -center`: the reference |
| `gromacs/md_protein.xtc`, `md_nowater.xtc` | the trajectory saved for part of the system, the topology for all of it |
| `gromacs/protein.gro` | the other way round: a protein-only structure, the full trajectory |
| `amber/complex.prmtop` + `md.nc` | Amber topology that carries PDB residue numbers and chains, NetCDF trajectory |
| `charmm/complex.psf` + `md.dcd` | CHARMM/NAMD topology with segment names, DCD trajectory |

```bash
bash examples/formats/run_all.sh          # twelve runs, tables only, a few minutes; exit code = runs that failed
python examples/formats/compare.py        # exit code 0 if every run is there and agrees with the reference
```

Every run must find the same 573 contacts between the same residues and atoms in the same frames, with
distances within 0.01 Å (xtc files round coordinates to 0.01 Å). `compare.py` prints one line per run:

```text
reference: 573 contacts in 11 frames, residues LYS 27 ... / ASP 39 ...
same      amber              573 contacts, missing 0, extra 0, largest distance difference 0.010 Å
same      charmm             573 contacts, missing 0, extra 0, largest distance difference 0.010 Å
...
same      waters             573 contacts, missing 0, extra 0, largest distance difference 0.010 Å, plus 211 water-mediated
```

The count of water-mediated contacts is 211 on our macOS machine and 212 on Linux in CI: one bridge sits
exactly on a cutoff. The 573 direct contacts are the same everywhere.

With MDAnalysis older than 2.11 the two `md.tpr` runs print `md.tpr is from a GROMACS newer than the installed
MDAnalysis reads: used md.gro next to it instead`. That is the documented behaviour, not a failure.

## 4. The Blender scene

```bash
conda create -n ppitango-blender python=3.11 && conda activate ppitango-blender
pip install ".[md,test,blender]"
pytest -q                                                   # now the Blender test runs too
ppitango md examples/simulations/output/md_1brs/top.pdb examples/simulations/output/md_1brs/traj.dcd \
    --chains A D --dt 0.1 --no-movie --no-viewer --blender -o out/md
python examples/blender_restyle.py out/md.blend out/restyled.png 200
```

`out/md.blend` should open in Blender 5.0 or newer, play with the space bar, and hold what the README's
Blender table lists. Nobody has yet clicked through the scene in the Blender window: a report from there is
especially welcome. Compare with [`demos/barnase_barstar_md.blend`](demos/barnase_barstar_md.blend).

## 5. Then try to break it

Ideas, roughly in order of how likely they are to find something:

- your own trajectories, from any MD program, any size, any box;
- other complexes from the PDB (`ppitango nma some.pdb`): antibodies with insertion codes (1HZH), more than
  two chains, a ligand or an ion as one side (`--select protein "resname ZN"`), nucleic acids (1ZAA);
- options the README lists: `--select`, `--renumber`, `--chains`, `--stride`, `--max-frames`, `--smooth`,
  `--waters`, `--max-residues`, `--seconds`, `--fps`, `--title`, `-o` into a folder that does not exist;
- wrong input on purpose: a chain that is not there, one chain only, a single frame, a trajectory that does
  not belong to the topology, an empty file. The README says such runs stop with a message that names the
  problem. A Python traceback is a finding;
- the Python API in the README's "From Python" section;
- Windows, Linux and macOS;
- every sentence under "What was tested, and known limits" in the README: is it true?

Already known and written down in the README's limits list, so not news: trajectories of a custom index
group; lone `.tpr` and plain `prmtop` files count residues straight through; a ligand must be a side on its
own; no 2D map or distance labels in the Blender scene; Blender needs Python 3.11.

## 6. Reporting

Open an issue at https://github.com/varunsharm16/ppi-tango/issues, one per finding, with:

1. the exact command, and where the input files come from (a path in this repository, or a PDB ID);
2. the complete output, traceback included;
3. what you expected, and the sentence of the README or of this file that made you expect it;
4. `ppitango --version`, `python --version`, the operating system, and `pip list` for MDAnalysis, numpy,
   pandas and pandappi.

Say how sure you are: "fails every time" and "looked odd in one frame" are both useful, but they are
different findings.
