"""Turn a ppi-tango MD run (and optionally an FEP result) into the data a live deck draws: two cartoons in one
frame, the Cα trace and the interface side chains along the trajectory, every contact with its occupancy,
and the closest polar distance of the residue pairs the talk zooms in on.

    python slides/tools/prep_interface.py examples/simulations/output/md_1brs/top.pdb \
        examples/simulations/output/md_1brs/traj.dcd --chains A D --dt 0.1 \
        --focus D:39 D:29 --dg examples/simulations/output/fep_1brs/dg.csv \
        slides/talks/barnase_barstar/data/interface.json

Writes <out>.json and <out>.js; the .js sets window.DATA["<out stem>"] (decks open from disk, where browsers
will not fetch JSON).

- Frame: Å, origin between the two faces of the interface, x from chain a to chain b (screen right), z along
  the interface's long direction (screen up), y into the screen: the turntable's starting view.
- Contacts are PandaPPI's, per frame, through ppitango.analyze. Occupancy is the fraction of frames a
  (type, residue, residue) contact is present.
- Drawn coordinates are a rolling mean over --smooth frames (thermal jitter makes a ribbon shiver); distances
  are measured on the raw frames.
- Secondary structure is DSSP through mdtraj on the first frame.

Needs ppi-tango[md] and mdtraj.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import ppitango

POLAR = ("N", "O")
CHARGE = {"ARG": 1, "LYS": 1, "HIS": 1, "ASP": -1, "GLU": -1}   # His counted as it sits against a carboxylate


def ca_frame(ca, i):
    """(tangent, bend, normal) at residue i; must match frameAt() in web/addons/cartoon.js."""
    i = min(max(i, 1), len(ca) - 2)
    t = ca[i + 1] - ca[i - 1]
    t /= np.linalg.norm(t)
    b = 2 * ca[i] - ca[i - 1] - ca[i + 1]
    b -= t * (t @ b)
    b /= np.linalg.norm(b)
    return t, b, np.cross(t, b)


def rolling(x, w):
    """Mean over a centred window of w frames along axis 0 (shorter at the ends)."""
    if w <= 1:
        return x
    out = np.empty_like(x)
    h = w // 2
    for f in range(len(x)):
        out[f] = x[max(0, f - h):f + h + 1].mean(0)
    return out


def label(resname, resnum):
    return f"{resname.title()}{resnum}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("top")
    ap.add_argument("traj")
    ap.add_argument("out", type=Path)
    ap.add_argument("--chains", nargs=2, required=True)
    ap.add_argument("--dt", type=float, default=1.0, help="time between trajectory frames, ns")
    ap.add_argument("--focus", nargs="*", default=[], help="residues the talk zooms in on, as CHAIN:NUMBER")
    ap.add_argument("--frames", type=int, default=100, help="frames kept for drawing")
    ap.add_argument("--smooth", type=int, default=5, help="rolling mean (frames) on drawn coordinates")
    ap.add_argument("--min-occupancy", type=float, default=0.2, help="contacts rarer than this are left out")
    ap.add_argument("--dg", help="FEP table from ppi-tango's example (window, lambdas, dG_complex, dG_free, dG, dG_err)")
    args = ap.parse_args()
    a, b = args.chains

    ens = ppitango.md(args.top, args.traj, chains=[a, b], dt=args.dt, max_frames=None)
    atoms, raw = ens.atoms.reset_index(drop=True), ens.coords
    n_frames = len(raw)
    result = ppitango.analyze(ens, progress=False)
    con = result.contacts
    con = con[~con.type.str.lower().str.contains("van der waals|vdw")]

    # ---- the frame ----
    key = lambda c, num: (c, int(num))
    res_of = {}
    for i, (c, num) in enumerate(zip(atoms.chain, atoms.resnum)):
        res_of.setdefault(key(c, num), []).append(i)
    num_of = lambda s: int(s.split()[1])
    face_a = np.concatenate([res_of[key(a, num_of(r))] for r in con.residue_a.unique()])
    face_b = np.concatenate([res_of[key(b, num_of(r))] for r in con.residue_b.unique()])
    ca_, cb_ = raw[0][face_a].mean(0), raw[0][face_b].mean(0)
    centre = (ca_ + cb_) / 2
    ex = (cb_ - ca_) / np.linalg.norm(cb_ - ca_)
    face = raw[0][np.concatenate([face_a, face_b])] - centre
    flat = face - np.outer(face @ ex, ex)
    ez = np.linalg.svd(flat, full_matrices=False)[2][0]
    ez -= ex * (ez @ ex)
    ez /= np.linalg.norm(ez)
    rot = np.stack([ex, np.cross(ez, ex), ez])
    raw = (raw - centre) @ rot.T
    drawn = rolling(raw, args.smooth)
    keep = np.unique(np.round(np.linspace(0, n_frames - 1, min(args.frames, n_frames))).astype(int))

    # ---- secondary structure ----
    import mdtraj
    t = mdtraj.load(args.top)
    prot = t.atom_slice(t.topology.select("protein"))
    dssp = mdtraj.compute_dssp(prot, simplified=True)[0]
    # mdtraj numbers chains by order; match residues by (order of appearance of chain, resSeq)
    order = []
    for line in open(args.top):
        if line.startswith(("ATOM", "HETATM")) and line[21] not in order:
            order.append(line[21])
    ss_of = {}
    for r, code in zip(prot.topology.residues, dssp):
        ss_of[(order[r.chain.index], r.resSeq)] = code if code in "HE" else "C"

    # ---- cartoons ----
    chains = {}
    for tag, c in (("a", a), ("b", b)):
        sel = atoms[atoms.chain == c]
        nums = list(dict.fromkeys(sel.resnum))
        at = lambda num, name: sel.index[(sel.resnum == num) & (sel.name == name)][0]
        ca_idx = [at(n, "CA") for n in nums]
        ca = raw[0][ca_idx]
        dirs = []
        for k, n in enumerate(nums):
            d = raw[0][at(n, "O")] - raw[0][at(n, "C")]
            d /= np.linalg.norm(d)
            tt, bb, nn = ca_frame(ca.copy(), k)
            dirs.append([round(float(d @ tt), 4), round(float(d @ bb), 4), round(float(d @ nn), 4)])
        names = [sel.resname[at(n, "CA")] for n in nums]
        chains[tag] = {
            "chain": c, "ca": np.round(ca, 2).tolist(), "dir": dirs, "breaks": [],
            "ss": "".join(ss_of.get((c, n), "C") for n in nums),
            "residues": [{"name": nm.title(), "num": int(n)} for nm, n in zip(names, nums)],
            "radius": round(float(np.linalg.norm(ca - ca.mean(0), axis=1).max()), 1),
            "centre": np.round(ca.mean(0), 2).tolist(),
            "trace": np.round(drawn[keep][:, ca_idx], 2).tolist(),
        }

    # ---- contacts ----
    groups = con.groupby(["type", "residue_a", "residue_b"])
    contacts = []
    for (typ, ra, rb), g in groups:
        present = np.zeros(n_frames, bool)
        present[g.frame.unique()] = True
        occ = present.mean()
        if occ < args.min_occupancy:
            continue
        contacts.append({
            "type": typ, "a": label(*ra.split()), "b": label(*rb.split()),
            "occupancy": round(float(occ), 3), "distance": round(float(g.distance.mean()), 2),
            "present": "".join("1" if p else "0" for p in present),
        })
    contacts.sort(key=lambda c: -c["occupancy"])

    # ---- residues on the map: every residue in a kept contact ----
    def node(c, lab):
        num = int("".join(ch for ch in lab if ch.isdigit()))
        idx = res_of[key(c, num)]
        ca_i = [i for i in idx if atoms.name[i] == "CA"][0]
        name = atoms.resname[ca_i]
        return {"id": lab, "name": name.title(), "num": num, "charge": CHARGE.get(name, 0),
                "at": np.round(raw[0][ca_i], 2).tolist()}
    nodes = {"a": [node(a, l) for l in dict.fromkeys(c["a"] for c in contacts)],
             "b": [node(b, l) for l in dict.fromkeys(c["b"] for c in contacts)]}

    # ---- residues to draw as sticks, and the pairs to measure ----
    focus = [(f.split(":")[0], int(f.split(":")[1])) for f in args.focus]
    polar_types = ("hydrogen bond", "salt bridge")
    pairs = []
    for c in contacts:
        if not any(p in c["type"].lower() for p in polar_types) or c["occupancy"] < 0.5:
            continue
        na, nb = (int("".join(ch for ch in c[k] if ch.isdigit())) for k in ("a", "b"))
        if ((a, na) in focus or (b, nb) in focus) and (c["a"], c["b"]) not in [(p["a"], p["b"]) for p in pairs]:
            pairs.append({"a": c["a"], "b": c["b"], "na": na, "nb": nb})
    shown = list(dict.fromkeys([(a, p["na"]) for p in pairs] + [(b, p["nb"]) for p in pairs] + focus))
    stick_idx, owner = [], []
    for k, (c, num) in enumerate(shown):
        for i in res_of[key(c, num)]:
            stick_idx.append(i)
            owner.append(k)
    local = {i: j for j, i in enumerate(stick_idx)}
    xyz0 = raw[0][stick_idx]
    bonds = [[i, j] for i in range(len(stick_idx)) for j in range(i + 1, len(stick_idx))
             if owner[i] == owner[j] and np.linalg.norm(xyz0[i] - xyz0[j]) < 1.9]
    sticks = {
        "residues": [{"chain": "a" if c == a else "b", "id": label(atoms.resname[res_of[key(c, num)][0]], num)}
                     for c, num in shown],
        "names": [atoms.name[i] for i in stick_idx], "elements": [atoms.element[i].title() for i in stick_idx],
        "residue": owner, "bonds": bonds, "xyz": np.round(drawn[keep][:, stick_idx], 2).tolist(),
    }
    for p in pairs:
        ia = [i for i in res_of[key(a, p.pop("na"))] if i in local and atoms.element[i] in POLAR]
        ib = [i for i in res_of[key(b, p.pop("nb"))] if i in local and atoms.element[i] in POLAR]
        d = np.linalg.norm(raw[:, ia][:, :, None] - raw[:, ib][:, None], axis=-1)        # frames, |ia|, |ib|
        best = d.reshape(n_frames, -1).argmin(1)
        p["i"] = [local[ia[k // len(ib)]] for k in best[keep]]
        p["j"] = [local[ib[k % len(ib)]] for k in best[keep]]
        p["d"] = np.round(d.reshape(n_frames, -1).min(1), 2).tolist()     # every raw frame, for the traces
        p["mean"], p["sd"] = round(float(np.mean(p["d"])), 2), round(float(np.std(p["d"])), 2)

    out = {
        "source": {"top": Path(args.top).name, "traj": Path(args.traj).name, "chains": [a, b]},
        "ns": round(float(ens.x[-1]), 2), "frames": n_frames, "kept": keep.tolist(),
        "time": np.round(np.asarray(ens.x)[keep], 2).tolist(),
        "chains": chains, "nodes": nodes, "contacts": contacts, "sticks": sticks, "pairs": pairs,
    }
    if args.dg:
        dg = pd.read_csv(args.dg)
        out["fep"] = {col: dg[col].round(3).tolist() for col in dg.columns}
    args.out.write_text(json.dumps(out, separators=(",", ":")))
    args.out.with_suffix(".js").write_text(
        f"(window.DATA = window.DATA || {{}})[{json.dumps(args.out.stem)}] = {json.dumps(out, separators=(',', ':'))};\n")
    always = sum(c["occupancy"] >= 0.95 for c in contacts)
    print(f"{n_frames} frames ({out['ns']} ns), {len(keep)} kept; {len(contacts)} contacts with occupancy >= "
          f"{args.min_occupancy}, {always} of them >= 0.95; {len(stick_idx)} stick atoms in "
          f"{len(shown)} residues; {len(pairs)} measured pairs: "
          + ", ".join(f"{p['a']}-{p['b']} {p['mean']}±{p['sd']} Å" for p in pairs))
    print(f"{args.out.with_suffix('.js')}: {args.out.with_suffix('.js').stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
