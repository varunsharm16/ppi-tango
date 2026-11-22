#!/bin/bash
# Rebuild everything in demos/. Run from the repository root.
# The MD and FEP inputs come from examples/simulations/simulate.py (see its docstring); the
# normal-mode demos need only the PDB files in examples/data.
set -e
D=examples/data
S=${SIM_DIR:-examples/simulations/output}   # where md_1brs/ and fep_1brs/ from simulate.py live

ppitango nma $D/1BRS.pdb --chains A D -o demos/barnase_barstar_nma \
    --title "Barnase–barstar (1BRS) · normal mode 1" --movie3d --gif
ppitango nma $D/1BJ1.pdb --chains H W -o demos/fab_vegf_nma \
    --title "Fab–VEGF (1BJ1) · normal mode 1" --movie3d
ppitango nma $D/6M0J.pdb --chains A E -o demos/ace2_rbd_nma \
    --title "ACE2–SARS-CoV-2 RBD (6M0J) · normal mode 1" --movie3d

ppitango md $S/md_1brs/top.pdb $S/md_1brs/traj.dcd --chains A D --dt 0.1 -o demos/barnase_barstar_md \
    --title "Barnase–barstar (1BRS) · 20 ns molecular dynamics" --seconds 16 --movie3d --gif --blender   # the .blend needs the blender extra (Python 3.11); skipped without
ppitango fep $S/fep_1brs/top.pdb $S/fep_1brs/window_*.dcd --chains A D --dg $S/fep_1brs/dg.csv --frames-per-window 10 \
    -o demos/barnase_barstar_fep --label "barstar Tyr29 side chain off" \
    --title "Barnase–barstar · FEP: barstar Tyr29 side chain switched off" --seconds 16 --movie3d
