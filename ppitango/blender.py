"""The Blender scene: the moving complex with its contacts and a camera path, saved as one .blend file to
restyle and render in Blender.

The models go in through Molecular Nodes, so its styles (cartoon, surface, sticks ...) are there to swap. The
contacts are one more node next to those styles: dashed lines between the two atoms of every contact, there
while the contact is made. Everything lives in the file itself; opening it needs neither ppi-tango nor the
trajectory.
"""
import os
import tempfile

import numpy as np
from pandappi.visualization.interaction_map import CHAIN_COLORS, _style

from . import movie
from .analysis import Scene
from .ensembles import _pdb_text
from .movie3d import TAN_HALF_FOV, _flight
from .viewer import _bridges, _canonical

SCALE = 0.01    # Molecular Nodes draws 1 Å as 0.01 Blender units
FPS = 24
WATER = '#E64B35'


def _linear(colour: str) -> list:
    """'#RRGGBB' as the linear RGBA Blender stores."""
    c = np.array([int(colour[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return [*np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4), 1.0]


def _contact_frames(scene: Scene, frames: list):
    """The contact lines of every frame, the same slots in the same order in each: end points (frames, slots,
    2, 3) in Å, occupancy (frames, slots), and per slot its colour and whether its second end is a bridging
    water. A water-mediated contact takes two slots, atom to water and water to atom."""
    ens = scene.ensemble
    coords, bridges = _canonical(scene), _bridges(scene)
    ends, occ, colours, water = [], [], [], []
    for k in (k for k in range(len(scene.contacts)) if scene.occ[k].max() >= 0.5):
        p = coords[frames, scene.atom_a[k, frames]]
        q = coords[frames, scene.atom_b[k, frames]]
        legs = [(p, q)]
        if k in bridges:
            w = bridges[k][frames]
            w = np.where(np.isnan(w), (p + q) / 2, w)       # no water where the contact is not made
            legs = [(p, w), (w, q)]
        for i, leg in enumerate(legs):
            ends.append(np.stack(leg, axis=1))
            occ.append(scene.occ[k, frames])
            colours.append(_linear(_style(scene.contacts[k]['type'])[0]))
            water.append(len(legs) == 2 and i == 0)
    return np.stack(ends, axis=1), np.stack(occ, axis=1), np.array(colours), np.array(water)


def _contacts_tree(bpy, n_frames: int, per_frame: int):
    """The node group that draws the contacts: picks the lines of the current frame out of the data object (all
    frames one after another, `per_frame` points each), blends towards the next frame, keeps the contacts that
    are made, and turns them into dashes."""
    tree = bpy.data.node_groups.new('ppi-tango contacts', 'GeometryNodeTree')
    for name, kind, value in (('Frame', 'Float', None), ('Contacts', 'Object', None), ('Shown above occupancy', 'Float', 0.5),
                              ('Radius (Å)', 'Float', 0.075), ('Dash (Å)', 'Float', 0.28), ('Gap (Å)', 'Float', 0.18),
                              ('Water radius (Å)', 'Float', 0.45), ('Material', 'Material', None)):
        socket = tree.interface.new_socket(name, in_out='INPUT', socket_type=f'NodeSocket{kind}')
        if value is not None:
            socket.default_value = value
    tree.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')

    def node(kind, inputs=None, **props):
        """A new node; inputs maps a socket (index or name) to a value or to the output socket feeding it."""
        n = tree.nodes.new(kind)
        for key, value in props.items():
            setattr(n, key, value)
        for key, value in (inputs or {}).items():
            if isinstance(value, bpy.types.NodeSocket):
                tree.links.new(value, n.inputs[key])
            else:
                n.inputs[key].default_value = value
        return n

    def math(operation, a, b=None):
        return node('ShaderNodeMath', {0: a} if b is None else {0: a, 1: b}, operation=operation).outputs[0]

    given = tree.nodes.new('NodeGroupInput').outputs
    scaled = lambda name: math('MULTIPLY', given[name], SCALE)
    data = node('GeometryNodeObjectInfo', {'Object': given['Contacts']}).outputs['Geometry']
    attribute = lambda name, kind: node('GeometryNodeInputNamedAttribute', {'Name': name},
                                        data_type=kind).outputs[0]
    # where we are: frame i, a fraction of the way to frame j
    u = node('ShaderNodeClamp', {'Value': given['Frame'], 'Min': 0.0, 'Max': n_frames - 1.0}).outputs[0]
    i = math('FLOOR', u)
    j = math('MINIMUM', math('ADD', i, 1.0), n_frames - 1.0)
    fraction = math('SUBTRACT', u, i)
    index = node('GeometryNodeInputIndex').outputs[0]

    def blended(value, kind):
        at = [node('GeometryNodeSampleIndex',
                   {'Geometry': data, 'Value': value, 'Index': math('ADD', index, math('MULTIPLY', frame, per_frame))},
                   data_type=kind, domain='POINT').outputs[0] for frame in (i, j)]
        mix = node('ShaderNodeMix', data_type='VECTOR' if kind == 'FLOAT_VECTOR' else 'FLOAT')
        tree.links.new(fraction, mix.inputs[0])
        a, b, out = (4, 5, 1) if kind == 'FLOAT_VECTOR' else (2, 3, 0)
        tree.links.new(at[0], mix.inputs[a])
        tree.links.new(at[1], mix.inputs[b])
        return mix.outputs[out]

    position = blended(node('GeometryNodeInputPosition').outputs[0], 'FLOAT_VECTOR')
    occupancy = blended(attribute('occupancy', 'FLOAT'), 'FLOAT')
    later = node('FunctionNodeCompare', {2: attribute('frame', 'INT'), 3: 0}, data_type='INT', operation='NOT_EQUAL')
    lines = node('GeometryNodeDeleteGeometry', {'Geometry': data, 'Selection': later.outputs[0]}, domain='POINT')
    lines = node('GeometryNodeSetPosition', {'Geometry': lines.outputs[0], 'Position': position})
    off = node('FunctionNodeCompare', {0: occupancy, 1: given['Shown above occupancy']}, data_type='FLOAT',
               operation='LESS_THAN')
    lines = node('GeometryNodeDeleteGeometry', {'Geometry': lines.outputs[0], 'Selection': off.outputs[0]},
                 domain='POINT').outputs[0]

    # dashes: a short cylinder at every step along each line, starting at its first atom
    dash = scaled('Dash (Å)')
    points = node('GeometryNodeCurveToPoints',
                  {'Curve': node('GeometryNodeMeshToCurve', {'Mesh': lines}).outputs[0],
                   'Length': math('ADD', dash, scaled('Gap (Å)'))}, mode='LENGTH')
    cylinder = node('GeometryNodeMeshCylinder', {'Vertices': 12, 'Radius': scaled('Radius (Å)'), 'Depth': dash})
    half = node('ShaderNodeCombineXYZ', {'Z': math('MULTIPLY', dash, 0.5)}).outputs[0]
    cylinder = node('GeometryNodeTransform', {'Geometry': cylinder.outputs['Mesh'], 'Translation': half})
    along = node('FunctionNodeAlignRotationToVector', {'Vector': points.outputs['Tangent']}, axis='Z')
    dashes = node('GeometryNodeInstanceOnPoints', {'Points': points.outputs['Points'],
                                                   'Instance': cylinder.outputs[0], 'Rotation': along.outputs[0]})
    # bridging waters: a sphere at the water end of the first leg
    waters = node('GeometryNodeMeshToPoints', {'Mesh': lines, 'Selection': attribute('water', 'BOOLEAN')})
    sphere = node('GeometryNodeMeshUVSphere', {'Segments': 24, 'Rings': 12, 'Radius': scaled('Water radius (Å)')})
    waters = node('GeometryNodeInstanceOnPoints', {'Points': waters.outputs[0], 'Instance': sphere.outputs['Mesh']})
    waters = node('GeometryNodeStoreNamedAttribute', {'Geometry': waters.outputs[0], 'Name': 'Color',
                                                      'Value': _linear(WATER)},
                  data_type='FLOAT_COLOR', domain='INSTANCE')

    joined = node('GeometryNodeJoinGeometry')
    for part in (dashes, waters):
        tree.links.new(part.outputs[0], joined.inputs[0])
    solid = node('GeometryNodeRealizeInstances', {'Geometry': joined.outputs[0]})
    solid = node('GeometryNodeSetShadeSmooth', {'Geometry': solid.outputs[0]})
    solid = node('GeometryNodeSetMaterial', {'Geometry': solid.outputs[0], 'Material': given['Material']})
    tree.links.new(solid.outputs[0], tree.nodes.new('NodeGroupOutput').inputs[0])
    return tree


def write(scene: Scene, output: str, seconds: float = 12.0) -> str:
    """Write the Blender scene of `scene` to `output` (.blend). Needs Blender as a Python module and Molecular
    Nodes: pip install "ppi-tango[blender]", on Python 3.11."""
    try:
        import bpy
        import molecularnodes as mn
        from molecularnodes.blender import coll
    except ImportError as e:
        raise RuntimeError(f'Blender or Molecular Nodes not found ({e}). Install them with: '
                           'pip install "ppi-tango[blender]"  (they need Python 3.11)')
    ens = scene.ensemble
    a, b = ens.chains
    atoms = ens.atoms
    n = max(int(seconds * FPS), 2)
    # a periodic motion ends on its first frame again, so that the last movie frame runs into the first
    frames = list(range(ens.n_frames)) + [0] * ens.periodic
    coords = _canonical(scene)

    bpy.ops.wm.read_factory_settings(use_empty=True)
    if not hasattr(bpy.types.Scene, 'MNSession'):
        mn.register()
    world = bpy.context.scene
    world.name = 'ppi-tango'
    world.frame_start, world.frame_end, world.render.fps = 1, n, FPS
    world.render.resolution_x, world.render.resolution_y = 1920, 1080
    world.view_settings.view_transform = 'Standard'         # colours as chosen, not toned down for film
    world.world = bpy.data.worlds.new('white')
    world.world.use_nodes = True                            # white to the camera, a dimmer grey as light
    sky = world.world.node_tree
    seen, lit = sky.nodes['Background'], sky.nodes.new('ShaderNodeBackground')
    seen.inputs['Color'].default_value, lit.inputs['Color'].default_value = (1, 1, 1, 1), (0.25, 0.25, 0.25, 1)
    mix = sky.nodes.new('ShaderNodeMixShader')
    for source, target in ((sky.nodes.new('ShaderNodeLightPath').outputs['Is Camera Ray'], mix.inputs[0]),
                           (lit.outputs[0], mix.inputs[1]), (seen.outputs[0], mix.inputs[2]),
                           (mix.outputs[0], sky.nodes['World Output'].inputs['Surface'])):
        sky.links.new(source, target)

    # the models, one per frame: Molecular Nodes keeps them in the file and moves the atoms between them
    with tempfile.TemporaryDirectory() as tmp:
        with open(f'{tmp}/frames.pdb', 'w') as f:
            f.write(''.join(f'MODEL {i + 1}\n' + _pdb_text(atoms, coords[frame]).replace('END\n', 'ENDMDL\n')
                            for i, frame in enumerate(frames)))
        mol = mn.Molecule.load(f'{tmp}/frames.pdb', name='complex', remove_solvent=False)
    if len(mol.object.data.vertices) != len(atoms):
        raise RuntimeError(f'Molecular Nodes read {len(mol.object.data.vertices)} atoms, ppi-tango has {len(atoms)}')
    shown = [k for k in range(len(scene.contacts)) if scene.occ[k].max() >= 0.5]
    at_interface = (set((a, scene.contacts[k]['residue_a']) for k in shown)
                    | set((b, scene.contacts[k]['residue_b']) for k in shown))
    interface = np.array([key in at_interface for key in zip(atoms.chain, ens.labels)])
    side_chains = interface & ~atoms.name.isin(['N', 'C', 'O']).values | (atoms.kind == 'other').values
    colour = np.array(mol.named_attribute('Color'))         # by element; carbons take the colour of their chain
    for chain, hue in zip((a, b), CHAIN_COLORS):
        colour[((atoms.chain == chain) & (atoms.element == 'C')).values] = _linear(hue)
    mol.store_named_attribute(colour, 'Color')
    mol.store_named_attribute(interface, 'interface')
    mol.store_named_attribute(side_chains, 'interface_side_chains')
    mol.add_style('cartoon', color=None)
    mol.add_style('ball_and_stick', color=None, selection='interface_side_chains')

    # one clock for everything: Molecular Nodes gives every style its own, keep the first
    tree = mol.tree
    movers = [n_ for n_ in tree.nodes if getattr(n_, 'node_tree', None) and n_.node_tree.name == 'Animate Frames']
    clocks = [n_ for n_ in tree.nodes if getattr(n_, 'node_tree', None) and n_.node_tree.name == 'Animate Value']
    for extra in movers[1:]:
        for link in list(extra.outputs[0].links):
            tree.links.new(movers[0].outputs[0], link.to_socket)
    for extra in movers[1:] + [c for c in clocks if not any(l.to_node == movers[0] for l in c.outputs[0].links)]:
        tree.nodes.remove(extra)
    clock = movers[0].inputs['Frame'].links[0].from_node
    clock.inputs['Frame Start'].default_value = 1
    clock.inputs['Frame End'].default_value = n + ens.periodic
    clock.inputs['Value Max'].default_value = len(frames) - 1
    clock.inputs['Clamped'].default_value = True

    # the contacts: every frame's lines in one hidden object, drawn by one more node next to the styles
    ends, occ, colours, water = _contact_frames(scene, frames)
    slots = ends.shape[1]
    mesh = bpy.data.meshes.new('contacts')
    mesh.from_pydata((ends * SCALE).reshape(-1, 3), np.arange(ends.size // 3).reshape(-1, 2), [])
    for name, kind, field, values in (
            ('occupancy', 'FLOAT', 'value', np.repeat(occ, 2)),
            ('frame', 'INT', 'value', np.repeat(np.arange(len(frames)), slots * 2)),
            ('water', 'BOOLEAN', 'value', np.tile(np.stack([np.zeros_like(water), water], 1).ravel(), len(frames))),
            ('Color', 'FLOAT_COLOR', 'color', np.tile(np.repeat(colours, 2, axis=0), (len(frames), 1)).ravel())):
        mesh.attributes.new(name, kind, 'POINT').data.foreach_set(field, values)
    data = bpy.data.objects.new('contacts', mesh)
    coll.data().objects.link(data)
    material = bpy.data.materials.new('contacts')
    material.use_nodes = True
    shader = material.node_tree
    shader.links.new(shader.nodes.new('ShaderNodeAttribute').outputs['Color'],
                     shader.nodes['Principled BSDF'].inputs['Base Color'])
    shader.nodes['Attribute'].attribute_name = 'Color'
    drawn = tree.nodes.new('GeometryNodeGroup')
    drawn.node_tree = _contacts_tree(bpy, len(frames), slots * 2)
    drawn.name = drawn.label = 'Contacts'
    drawn.inputs['Contacts'].default_value = data
    drawn.inputs['Material'].default_value = material
    tree.links.new(clock.outputs[0], drawn.inputs['Frame'])
    join = next(n_ for n_ in tree.nodes if n_.bl_idname == 'GeometryNodeJoinGeometry')
    tree.links.new(drawn.outputs[0], join.inputs[0])

    # the stage: chain a left, chain b right, the interface upright; the camera on a pivot that flies in to
    # the interface, rocks to either side of it and flies out again
    stage = bpy.data.objects.new('stage', None)
    stage.rotation_euler = (np.pi / 2, 0, 0)                # the upright axis of the models is Blender's z
    pivot = bpy.data.objects.new('camera pivot', None)
    camera = bpy.data.objects.new('camera', bpy.data.cameras.new('camera'))
    sun = bpy.data.objects.new('light', bpy.data.lights.new('light', 'SUN'))
    sun.data.energy = 2.5
    for obj, parent in ((stage, None), (pivot, stage), (camera, pivot), (sun, camera)):
        world.collection.objects.link(obj)
        obj.parent = parent
    mol.object.parent = stage
    world.camera = camera
    camera.data.sensor_fit = 'VERTICAL'
    camera.data.angle_y = 2 * np.arctan(TAN_HALF_FOV)
    reach = float(np.linalg.norm(coords[0], axis=1).max()) + 8
    face = coords[0][np.unique(np.concatenate([scene.atom_a[shown].ravel(), scene.atom_b[shown].ravel()]))]
    near = max(float(np.percentile(np.linalg.norm(face, axis=1), 90)) * 1.25, 14.0)
    far = reach * 0.8
    camera.data.clip_start, camera.data.clip_end = SCALE, 4 * far / TAN_HALF_FOV * SCALE
    for frame in sorted({*range(1, n, FPS // 4), n}):       # a key every quarter second, to edit by hand
        theta, tilt, zoom = _flight((frame - 1) / (n - 1))
        pivot.rotation_euler = (-tilt, theta, 0)
        camera.location = (0, 0, (far + (near - far) * zoom) / TAN_HALF_FOV * SCALE)
        pivot.keyframe_insert('rotation_euler', frame=frame)
        camera.keyframe_insert('location', frame=frame)
    world.frame_set(1)

    # the scene and what it uses, nothing else: a plain save also stores where Blender and Molecular Nodes
    # are installed on this computer
    bpy.data.libraries.write(os.path.abspath(output), {world}, compress=True)
    return output
