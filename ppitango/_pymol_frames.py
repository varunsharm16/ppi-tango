"""Run inside PyMOL (pymol -cq _pymol_frames.py -- job.json part parts): ray-trace the movie frames.

The job file holds a multi-state PDB, and for every frame the camera and the contacts to draw.
"""
import json
import sys

from pymol import cgo, cmd

job = json.load(open(sys.argv[1]))
part, parts = int(sys.argv[2]), int(sys.argv[3])
a, b = job['chains']
rgb = lambda h: [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)]

cmd.load(job['pdb'], 'm')
cmd.dss('m', state=1)
cmd.hide('everything')
cmd.bg_color('white')
for chain, name in ((a, 'tango_a'), (b, 'tango_b')):
    cmd.set_color(name, rgb(job['colors'][chain]))
    cmd.show('cartoon', f'm and chain {chain}')
    side = f"m and chain {chain} and resi {'+'.join(map(str, job['interface'][chain]))} and not name N+C+O"
    cmd.show('sticks', side)
    cmd.color(name, f'm and chain {chain}')
cmd.color('atomic', 'm and not elem C')
for key, value in dict(stick_radius=0.2, ray_trace_mode=1, ray_trace_gain=0.05,
                       antialias=1, ray_shadows=0, ambient=0.55, specular=0.15, ray_opaque_background=1,
                       cartoon_fancy_helices=1,
                       depth_cue=1, fog_start=0.55, ray_trace_fog=1, field_of_view=14).items():
    cmd.set(key, value)

DASH, GAP, RADIUS = 0.28, 0.18, 0.075
for i, frame in enumerate(job['frames']):
    if i % parts != part:
        continue
    obj = []
    for x1, y1, z1, x2, y2, z2, color in frame['contacts']:
        r, g, bl = rgb(color)
        length = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
        t = 0.0
        while t < length:
            u, v = t / length, min(t + DASH, length) / length
            obj += [cgo.CYLINDER, x1 + (x2 - x1) * u, y1 + (y2 - y1) * u, z1 + (z2 - z1) * u,
                    x1 + (x2 - x1) * v, y1 + (y2 - y1) * v, z1 + (z2 - z1) * v, RADIUS, r, g, bl, r, g, bl]
            t += DASH + GAP
    for x, y, z in frame.get('waters', ()):          # bridging waters: red spheres
        obj += [cgo.COLOR, 0.90, 0.29, 0.21, cgo.SPHERE, x, y, z, 0.45]
    if obj:
        cmd.load_cgo(obj, 'contacts', state=i + 1)

for i, frame in enumerate(job['frames']):
    if i % parts != part:
        continue
    cmd.set('state', i + 1)
    cmd.set_view(frame['view'])
    cmd.png(f"{job['out']}/f{i:05d}.png", width=job['width'], height=job['height'], ray=1, quiet=1)
