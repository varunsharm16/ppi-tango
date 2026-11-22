#!/bin/bash
# One system (barnase-barstar in water, 11 frames over 20 ps), saved the ways MD programs save it.
# Every line must find the same contacts. Run from the repository root:
#     bash examples/formats/run_all.sh        then        python examples/formats/compare.py
# Tables only (no movies), about a minute per line on a laptop. Results go to formats_out/.
G=examples/formats/gromacs; A=examples/formats/amber; C=examples/formats/charmm; O=formats_out
rm -rf $O   # tables of an earlier run must not be taken for this one's
failed=0
run() { name=$1; shift; echo "== $name: ppitango md $*"
        ppitango md "$@" --no-movie --no-viewer -o $O/$name || { echo "FAILED: $name"; failed=$((failed + 1)); }; }

# GROMACS: the whole system (37044 atoms: water, ions, dodecahedron box). md.xtc is as mdrun wrote it,
# both chains cut by the box in every frame; md_whole.xtc is the same after trjconv -pbc mol -center.
run reference        $G/md.gro    $G/md_whole.xtc
run gro_raw          $G/md.gro    $G/md.xtc
run tpr_raw          $G/md.tpr    $G/md.xtc
run top_raw          $G/topol.top $G/md.xtc
# the trajectory saved for part of the system, against the topology of all of it
run gro_protein_xtc  $G/md.gro    $G/md_protein.xtc
run tpr_protein_xtc  $G/md.tpr    $G/md_protein.xtc
run top_nowater_xtc  $G/topol.top $G/md_nowater.xtc
# a protein-only structure file, against the full and the protein-only trajectory
run protein_gro_full $G/protein.gro $G/md.xtc
run protein_gro      $G/protein.gro $G/md_protein.xtc
# Amber and CHARMM/NAMD: solute only, still cut by the box. The prmtop carries the PDB residue numbers and
# chains (parmed addPDB), the psf its resids and segment names.
run amber            $A/complex.prmtop $A/md.nc
run charmm           $C/complex.psf    $C/md.dcd
# water-mediated contacts on top (needs the water: full topology and full trajectory)
run waters           $G/md.gro    $G/md.xtc --waters

[ $failed -eq 0 ] || echo "$failed of 12 runs failed"
exit $failed
