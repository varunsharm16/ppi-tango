"""Compare the tables run_all.sh wrote: every run must hold the contacts of the reference run.

    python examples/formats/compare.py [formats_out]
"""
import os
import sys

import pandas as pd

out = sys.argv[1] if len(sys.argv) > 1 else 'formats_out'
key = ['frame', 'type', 'residue_a', 'atom_a', 'residue_b', 'atom_b']
load = lambda name: pd.read_csv(f'{out}/{name}_contacts.csv')
if not os.path.exists(f'{out}/reference_contacts.csv'):
    sys.exit(f'{out}/reference_contacts.csv is missing: run bash examples/formats/run_all.sh first')
ref = load('reference')
want = set(map(tuple, ref[key].values))
print(f'reference: {len(ref)} contacts in {ref.frame.nunique()} frames, residues '
      f'{ref.residue_a.iloc[0]} ... / {ref.residue_b.iloc[0]} ...')
RUNS = ('amber', 'charmm', 'gro_protein_xtc', 'gro_raw', 'protein_gro', 'protein_gro_full', 'top_nowater_xtc',
        'top_raw', 'tpr_protein_xtc', 'tpr_raw', 'waters')       # every line of run_all.sh but the reference
bad = 0
for name in RUNS:
    if not os.path.exists(f'{out}/{name}_contacts.csv'):       # a run that failed counts, it is not skipped
        print(f'MISSING   {name:18s} no table: the run failed or was not made')
        bad += 1
        continue
    df = load(name)
    if name == 'waters':                                     # the same contacts, plus the water-mediated ones
        extra = df[df.type.str.contains('ater')]
        df = df[~df.type.str.contains('ater')]
    got = set(map(tuple, df[key].values))
    order = lambda d: d.fillna('').sort_values(key + ['distance']).distance.values   # atom_b is empty for ring contacts
    worst = abs(order(ref) - order(df)).max() if len(df) == len(ref) else float('nan')
    same = got == want
    bad += not same
    note = f', plus {len(extra)} water-mediated' if name == 'waters' else ''
    print(f"{'same     ' if same else 'DIFFERENT'} {name:18s} {len(df)} contacts, missing {len(want - got)}, "
          f"extra {len(got - want)}, largest distance difference {worst:.3f} Å{note}")
sys.exit(bad)
