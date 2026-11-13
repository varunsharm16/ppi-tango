"""Simulations behind the ppi-tango demos (OpenMM, explicit solvent, Amber14SB/TIP3P).

    python simulate.py md  1BRS.pdb A D --ns 20 --out md_1brs
    python simulate.py fep 1BRS.pdb A D --mutate D:29 --leg complex --out fep_1brs
    python simulate.py fep 1BRS.pdb A D --mutate D:29 --leg free    --out fep_1brs
    python simulate.py fep-analyze fep_1brs

`fep` annihilates one side chain beyond CB (a pseudo-alanine mutation) along 12 lambda windows,
once in the complex and once in the isolated partner; `fep-analyze` runs MBAR on both legs and
writes dg.csv (cumulative free energy per window) for `ppitango fep`.
Needs openmm, pdbfixer; fep also needs openmmtools and pymbar.
"""
import argparse
import json
import os

import numpy as np
import openmm as mm
from openmm import app, unit

KT = 0.593  # kcal/mol at 298 K
# (lambda_electrostatics, lambda_sterics): charges off first, then soft-core Lennard-Jones
SCHEDULE = [(1, 1), (.75, 1), (.5, 1), (.25, 1), (0, 1),
            (0, .85), (0, .7), (0, .55), (0, .4), (0, .25), (0, .1), (0, 0)]


def build(pdb, chains, hmr=True):
    """Solvated system for the given chains. Returns (topology, positions, system, protein atom indices)."""
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=pdb)
    fixer.removeChains([i for i, c in enumerate(fixer.topology.chains()) if c.id not in chains])
    # the same chain id can be reused for waters/ligands: keep the first (polymer) chain of each id
    seen, drop = set(), []
    for i, c in enumerate(fixer.topology.chains()):
        if c.id in seen:
            drop.append(i)
        seen.add(c.id)
    fixer.removeChains(drop)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # do not model unresolved loops/termini
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    ff = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    model = app.Modeller(fixer.topology, fixer.positions)
    model.addSolvent(ff, padding=1.0 * unit.nanometer, ionicStrength=0.15 * unit.molar)
    system = ff.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds, hydrogenMass=(3 if hmr else 1) * unit.amu)
    system.addForce(mm.MonteCarloBarostat(1 * unit.bar, 298 * unit.kelvin))
    protein = [a.index for a in model.topology.atoms() if a.residue.chain.id in chains
               and a.residue.name not in ('HOH', 'NA', 'CL')]
    return model.topology, model.positions, system, protein


def simulation(topology, system, positions, dt_fs):
    integrator = mm.LangevinMiddleIntegrator(298 * unit.kelvin, 1 / unit.picosecond, dt_fs * unit.femtosecond)
    sim = app.Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    return sim


def write_top(topology, positions, protein, path):
    """Protein-only PDB matching the atom order of the saved trajectories."""
    model = app.Modeller(topology, positions)
    keep = set(protein)
    model.delete([a for a in model.topology.atoms() if a.index not in keep])
    with open(path, 'w') as f:
        app.PDBFile.writeFile(model.topology, model.positions, f, keepIds=True)


def protein_dcd(sim, protein, path, every):
    import mdtraj
    sim.reporters.append(mdtraj.reporters.DCDReporter(path, every, atomSubset=protein))


def run_md(args):
    os.makedirs(args.out, exist_ok=True)
    topology, positions, system, protein = build(args.pdb, args.chains)
    sim = simulation(topology, system, positions, 4)
    sim.minimizeEnergy()
    write_top(topology, sim.context.getState(getPositions=True).getPositions(), protein, f'{args.out}/top.pdb')
    sim.context.setVelocitiesToTemperature(298 * unit.kelvin)
    sim.step(50_000)  # 200 ps equilibration
    steps = int(args.ns * 250_000)
    protein_dcd(sim, protein, f'{args.out}/traj.dcd', steps // args.frames)
    sim.reporters.append(app.StateDataReporter(f'{args.out}/md.log', 25_000, time=True, potentialEnergy=True,
                                               temperature=True, speed=True))
    sim.step(steps)
    json.dump({'ns': args.ns, 'frames': args.frames, 'ns_per_frame': args.ns / args.frames},
              open(f'{args.out}/meta.json', 'w'))


def run_fep(args):
    from openmmtools import alchemy
    os.makedirs(args.out, exist_ok=True)
    chain, resnum = args.mutate.split(':')
    chains = args.chains if args.leg == 'complex' else [chain]
    topology, positions, system, protein = build(args.pdb, chains, hmr=False)
    backbone = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB2', 'HB3', 'H1', 'H2', 'H3', 'OXT'}
    side = [a.index for a in topology.atoms() if a.residue.chain.id == chain
            and a.residue.id == resnum and a.name not in backbone]
    resname = next(r.name for r in topology.residues() if r.chain.id == chain and r.id == resnum)
    print(f'{args.leg}: annihilating {len(side)} side-chain atoms of {chain}:{resname}{resnum}', flush=True)

    region = alchemy.AlchemicalRegion(alchemical_atoms=side, annihilate_electrostatics=True,
                                      annihilate_sterics=False)
    alch_system = alchemy.AbsoluteAlchemicalFactory().create_alchemical_system(system, region)
    state = alchemy.AlchemicalState.from_system(alch_system)
    sim = simulation(topology, alch_system, positions, 2)
    sim.minimizeEnergy()
    if args.leg == 'complex':
        write_top(topology, sim.context.getState(getPositions=True).getPositions(), protein, f'{args.out}/top.pdb')
    sim.context.setVelocitiesToTemperature(298 * unit.kelvin)

    def set_lambda(k):
        state.lambda_electrostatics, state.lambda_sterics = SCHEDULE[k]
        state.apply_to_context(sim.context)

    n_samples, every = int(args.ns * 1000), 500  # one energy sample per ps
    K = len(SCHEDULE)
    u = np.zeros((K, K, n_samples))  # u[k, l, n]: sample n of window k evaluated in state l (kT)
    kT = unit.MOLAR_GAS_CONSTANT_R * 298 * unit.kelvin
    for k in range(K):
        set_lambda(k)
        sim.step(int(os.environ.get("TANGO_EQUIL", 50_000)))  # 100 ps equilibration in each window
        sim.reporters.clear()
        if args.leg == 'complex':
            protein_dcd(sim, protein, f'{args.out}/window_{k:02d}.dcd', every * n_samples // args.frames)
        for n in range(n_samples):
            sim.step(every)
            for l in range(K):
                set_lambda(l)
                st = sim.context.getState(getEnergy=True)
                volume = st.getPeriodicBoxVolume() * 1 * unit.bar * unit.AVOGADRO_CONSTANT_NA
                u[k, l, n] = (st.getPotentialEnergy() + volume) / kT
            set_lambda(k)
        print(f'window {k} {SCHEDULE[k]} done', flush=True)
        np.save(f'{args.out}/u_{args.leg}.npy', u)
    json.dump({'mutate': args.mutate, 'resname': resname, 'ns_per_window': args.ns, 'schedule': SCHEDULE},
              open(f'{args.out}/meta.json', 'w'))


def mbar(u):
    """Cumulative free energy (kcal/mol) and its uncertainty along the windows."""
    from pymbar import MBAR, timeseries
    K = u.shape[0]
    cols, N = [], []
    for k in range(K):
        neighbour = k + 1 if k + 1 < K else k - 1
        g = timeseries.statistical_inefficiency(u[k, neighbour] - u[k, k])
        idx = timeseries.subsample_correlated_data(u[k, k], g=g)
        cols.append(u[k][:, idx])
        N.append(len(idx))
    res = MBAR(np.concatenate(cols, axis=1), N).compute_free_energy_differences()
    return res['Delta_f'][0] * KT, res['dDelta_f'][0] * KT


def analyze_fep(args):
    meta = json.load(open(f'{args.out}/meta.json'))
    (gc, ec), (gf, ef) = mbar(np.load(f'{args.out}/u_complex.npy')), mbar(np.load(f'{args.out}/u_free.npy'))
    with open(f'{args.out}/dg.csv', 'w') as f:
        f.write('window,lambda_electrostatics,lambda_sterics,dG_complex,dG_free,dG,dG_err\n')
        for k, (le, ls) in enumerate(meta['schedule']):
            f.write(f'{k},{le},{ls},{gc[k]:.3f},{gf[k]:.3f},{gc[k] - gf[k]:.3f},{np.hypot(ec[k], ef[k]):.3f}\n')
    print(f"ddG_bind({meta['mutate']} {meta['resname']} side chain off) = "
          f"{gc[-1] - gf[-1]:.2f} +- {np.hypot(ec[-1], ef[-1]):.2f} kcal/mol "
          f"(complex {gc[-1]:.2f}, free {gf[-1]:.2f})")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    for name in ('md', 'fep'):
        s = sub.add_parser(name)
        s.add_argument('pdb')
        s.add_argument('chains', nargs='+')
        s.add_argument('--out', required=True)
        s.add_argument('--ns', type=float, default=20 if name == 'md' else 1,
                       help='production length (per window for fep)')
        s.add_argument('--frames', type=int, default=200 if name == 'md' else 10,
                       help='saved frames (per window for fep)')
        if name == 'fep':
            s.add_argument('--mutate', required=True, help='chain:residue number, e.g. D:29')
            s.add_argument('--leg', choices=('complex', 'free'), required=True)
    sub.add_parser('fep-analyze').add_argument('out')
    args = p.parse_args()
    {'md': run_md, 'fep': run_fep, 'fep-analyze': analyze_fep}[args.cmd](args)
