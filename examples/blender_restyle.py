"""Restyle a ppi-tango Blender scene from a script instead of by hand, and render it.

    python examples/blender_restyle.py demos/barnase_barstar_md.blend restyled.mp4

Everything done here can be done in Blender itself: open the .blend, select the object "complex" and look
at its Geometry Nodes. This script wraps the cartoon in a see-through surface, makes the contacts thicker and
glowing, and turns the backdrop dark. Needs the same install as `ppitango --blender`.
"""
import os
import sys

import bpy
from molecularnodes.nodes import nodes

blend, output = sys.argv[1], os.path.abspath(sys.argv[2])
bpy.ops.wm.open_mainfile(filepath=blend)
scene = bpy.context.scene
tree = bpy.data.objects['complex'].modifiers[0].node_group

# a see-through surface around the cartoon: one more Molecular Nodes style, fed by the same moving atoms
cartoon = next(n for n in tree.nodes if n.bl_idname == 'GeometryNodeGroup' and n.node_tree.name == 'Style Cartoon')
glass = bpy.data.materials.new('glass')
glass.use_nodes = True
shader = glass.node_tree
colour = shader.nodes.new('ShaderNodeAttribute')
colour.attribute_name = 'Color'
shader.links.new(colour.outputs['Color'], shader.nodes['Principled BSDF'].inputs['Base Color'])
shader.nodes['Principled BSDF'].inputs['Alpha'].default_value = 0.22
glass.surface_render_method = 'BLENDED'
glass.use_backface_culling = True
surface = tree.nodes.new('GeometryNodeGroup')
surface.node_tree = nodes.append('Style Surface')
surface.inputs['Material'].default_value = glass
tree.links.new(cartoon.inputs[0].links[0].from_socket, surface.inputs[0])
tree.links.new(surface.outputs[0], cartoon.outputs[0].links[0].to_socket)

# the contacts: thicker, and glowing in their own colours
tree.nodes['Contacts'].inputs['Radius (Å)'].default_value = 0.12
shader = bpy.data.materials['contacts'].node_tree
shader.links.new(shader.nodes['Attribute'].outputs['Color'], shader.nodes['Principled BSDF'].inputs['Emission Color'])
shader.nodes['Principled BSDF'].inputs['Emission Strength'].default_value = 3

# a dark backdrop (the world holds two colours: what the camera sees, and what lights the models)
for node in scene.world.node_tree.nodes:
    if node.bl_idname == 'ShaderNodeBackground':
        seen = node.inputs['Color'].default_value[0] == 1
        node.inputs['Color'].default_value = (0.012, 0.016, 0.035, 1) if seen else (0.35, 0.35, 0.4, 1)

if output.endswith('.png'):                                  # one picture, of the frame given third
    scene.frame_set(int(sys.argv[3]) if len(sys.argv) > 3 else scene.frame_end // 2)
    scene.render.filepath = output
    bpy.ops.render.render(write_still=True)
else:                                                        # the whole movie
    scene.render.image_settings.media_type = 'VIDEO'
    scene.render.ffmpeg.format, scene.render.ffmpeg.codec = 'MPEG4', 'H264'
    scene.render.filepath = output
    bpy.ops.render.render(animation=True)
