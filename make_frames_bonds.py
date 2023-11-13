import bpy
import bmesh
import numpy as np
from os import scandir
import sys
sys.path.insert(0, r"C:\Users\jacks\Desktop\Blender\python")
import import_dump


# define input ============================================================================================
dump_location = r"test_dumps"  # folder if individual. file if composite
composite = False
coloring_field = 'v_vdia'.upper()  # name of heading to color by. Can change later in geonodes. must be uppercase
bond_id1_ix = 0  # index of atom id 1 in bond dumps
bond_id2_ix = bond_id1_ix + 1  # index of atom id 2 in bond dumps
bond_dx_ix = 2  # index of x component of bond position vector in bond dumps. x, y must folow
bond_dist_ix = 5  # index of bond distance in bond dumps
sortkey = import_dump.dumpnum  # how should the dumps be sorted
# define input ============================================================================================

print("loading dump file")
if composite:
    raise Exception("no composites please")
    dump_data_list = import_dump.lammps_composite(dump_location)
    fields, atoms, N, time = dump_data_list[0]
else:
    dump_data_list = [f.path for f in scandir(dump_location) if "dump" in f.name and "atom" in f.name]
    bond_data_list = [f.path for f in scandir(dump_location) if "dump" in f.name and "bond" in f.name]
    dump_data_list.sort(key=import_dump.dumpnum)
    bond_data_list.sort(key=sortkey)
    fields, atoms, N, time = import_dump.lammps_single(dump_data_list[0])
    bfields, bonds, bN, _ = import_dump.lammps_bond_single(bond_data_list[0], cutoff=2.0, cutoff_index=bond_dist_ix)


bpy.context.scene.frame_end = len(dump_data_list)
# define meshes
print("creating meshes")

mesh = bpy.data.meshes.new("mesh")
atom_obj = bpy.data.objects.new('atoms', mesh)
bpy.context.scene.collection.children[0].objects.link(atom_obj)

bond_mesh = bpy.data.meshes.new("mesh")
bond_obj = bpy.data.objects.new('bonds', bond_mesh)
bpy.context.scene.collection.children[0].objects.link(bond_obj)

# create mesh geometry data
print("defining geometry data")
if "ux" in fields:
    x_ix = fields.index("ux")
    y_ix = fields.index("uy")
    z_ix = fields.index("uz")
else:
    x_ix = fields.index("x")
    y_ix = fields.index("y")
    z_ix = fields.index("z")
rs = np.array([[atom[x_ix], atom[y_ix], atom[z_ix]] for atom in atoms])
mesh.from_pydata(rs, [], [])

if np.all(np.arange(1, N + 1) == np.array([a[0] for a in atoms])):
        midpoints = [np.array(atoms[int(np.minimum(b[bond_id1_ix],b[bond_id2_ix])-1)][x_ix:z_ix+1]) - np.array(b[bond_dx_ix:bond_dx_ix+3])/2 for b in bonds]
else:
      atom_ids = [a[0] for a in atoms]
      midpoints = [np.array(atoms[atom_ids.index(np.minimum(b[bond_id1_ix],b[bond_id2_ix]))][x_ix:z_ix+1]) - np.array(b[bond_dx_ix:bond_dx_ix+3])/2 for b in bonds]
#       raise Exception("implement bond endpoint finding for non-sequential atom numbering")
rotations = [import_dump.vec2rot(b[bond_dx_ix:bond_dx_ix+3]) for b in bonds]
scales = [b[bond_dist_ix] for b in bonds]
bond_mesh.from_pydata(midpoints,[], [])

bond_obj.data.attributes.new(name='rotation', type='FLOAT_VECTOR', domain='POINT')
rotations = np.array(rotations).reshape(-1)
bond_obj.data.attributes['rotation'].data.foreach_set('vector', rotations)

bond_obj.data.attributes.new(name='scale', type='FLOAT', domain='POINT')
bond_obj.data.attributes['scale'].data.foreach_set('value', scales)

print("defining coloring data")
# define coloring attribute data
for i, field in enumerate(fields):
        atom_obj.data.attributes.new(name=field.upper(), type='FLOAT', domain='POINT')
        atom_obj.data.attributes[field.upper()].data.foreach_set('value', [atom[i] for atom in atoms])
for i, field in enumerate(bfields):
        bond_obj.data.attributes.new(name=field.upper(), type='FLOAT', domain='POINT')
        bond_obj.data.attributes[field.upper()].data.foreach_set('value', [bond[i] for bond in bonds])


print("defining materials")

# make valued materials
atom_value_material = bpy.data.materials.new(name="atom_val_mat")
atom_value_material.use_nodes = True
bond_value_material = bpy.data.materials.new(name="bond_val_mat")
bond_value_material.use_nodes = True
# make chemical materials
atom_1_material = bpy.data.materials.new(name="atom_1_mat")
atom_1_material.use_nodes = True
atom_2_material = bpy.data.materials.new(name="atom_2_mat")
atom_2_material.use_nodes = True

print("defining atom modifier")
# geometry nodes
    # atoms
bpy.context.view_layer.objects.active = atom_obj
bpy.ops.object.modifier_add(type='NODES')
bpy.ops.node.new_geometry_node_group_assign()
ng = atom_obj.modifiers[-1].node_group
nodes = ng.nodes

input = nodes.get('Group Input')
input.location = (-1500, 0)
output = nodes.get('Group Output')
output.location = (1000, 0)
        # attributes
print("defining atom geonodes attributes")
new_attribute = ng.inputs.new('NodeSocketFloat', "TYPE")
atom_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'TYPE'
atom_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
new_attribute = ng.inputs.new('NodeSocketFloat', coloring_field)
atom_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = coloring_field
atom_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
new_attribute = ng.outputs.new('NodeSocketColor','AtomColor')
atom_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'atom_color'
atom_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
        # nodes
print("defining atom nodes")        
sphere_node = nodes.new('GeometryNodeMeshUVSphere')
sphere_node.location = (-1000, -600)
sphere_node.inputs['Radius'].default_value=0.40
set_material_node = nodes.new('GeometryNodeSetMaterial')
set_material_node.location = (-600, -600)
set_material_node.inputs['Material'].default_value = atom_value_material
smooth_node = nodes.new('GeometryNodeSetShadeSmooth')
smooth_node.location = (-800,-600)
instance_node = nodes.new('GeometryNodeInstanceOnPoints')
instance_node.location = (-400, -500)

sphere_node_1 = nodes.new('GeometryNodeMeshUVSphere')
sphere_node_1.location = (-1000, 200)
sphere_node_1.inputs['Radius'].default_value=0.2
set_material_node_1 = nodes.new('GeometryNodeSetMaterial')
set_material_node_1.location = (-600, 200)
set_material_node_1.inputs['Material'].default_value = atom_1_material
smooth_node_1 = nodes.new('GeometryNodeSetShadeSmooth')
smooth_node_1.location = (-800, 200)
instance_node_1 = nodes.new('GeometryNodeInstanceOnPoints')
instance_node_1.location = (-400, 300)

sphere_node_2 = nodes.new('GeometryNodeMeshUVSphere')
sphere_node_2.location = (-1000, -200)
sphere_node_2.inputs['Radius'].default_value=0.6
set_material_node_2 = nodes.new('GeometryNodeSetMaterial')
set_material_node_2.location = (-600, -200)
set_material_node_2.inputs['Material'].default_value = atom_2_material
smooth_node_2 = nodes.new('GeometryNodeSetShadeSmooth')
smooth_node_2.location = (-800, -200)
instance_node_2 = nodes.new('GeometryNodeInstanceOnPoints')
instance_node_2.location = (-400, -100)

minmax_node = nodes.new('GeometryNodeAttributeStatistic')
minmax_node.location = (0, -400)
minmax_node.domain = 'INSTANCE' 
map_range_node = nodes.new('ShaderNodeMapRange')
map_range_node.location = (200, -400)
color_node = nodes.new('ShaderNodeValToRGB')
color_node.location = (400, -400)

ID_compare_node_1 = nodes.new('FunctionNodeCompare')
ID_compare_node_1.location = (-1000, 400)
ID_compare_node_1.inputs[1].default_value=2
ID_compare_node_1.operation = 'EQUAL'
ID_compare_node_1.data_type = 'FLOAT'
separate_geometry_node_1 = nodes.new('GeometryNodeSeparateGeometry')
separate_geometry_node_1.location = (-800, 400)
separate_geometry_node_1.domain = 'POINT'

ID_compare_node_2 = nodes.new('FunctionNodeCompare')
ID_compare_node_2.location = (-1000, 000)
ID_compare_node_2.inputs[1].default_value=3
ID_compare_node_2.operation = 'EQUAL'
ID_compare_node_2.data_type = 'FLOAT'
separate_geometry_node_2 = nodes.new('GeometryNodeSeparateGeometry')
separate_geometry_node_2.location = (-800, 000)
separate_geometry_node_2.domain = 'POINT'

join_geometry_node = nodes.new('GeometryNodeJoinGeometry')
join_geometry_node.location = (600, 0)
realize_instances_node = nodes.new('GeometryNodeRealizeInstances')
realize_instances_node.location = (800, 0)
        #links: geometry
print("defining atom links")
ng.links.new(input.outputs['Geometry'], separate_geometry_node_1.inputs['Geometry'])
ng.links.new(separate_geometry_node_1.outputs['Selection'], instance_node_1.inputs['Points'])
ng.links.new(separate_geometry_node_1.outputs['Inverted'], separate_geometry_node_2.inputs['Geometry'])

ng.links.new(separate_geometry_node_2.outputs['Selection'], instance_node_2.inputs['Points'])
ng.links.new(separate_geometry_node_2.outputs['Inverted'], instance_node.inputs['Points'])

ng.links.new(sphere_node.outputs['Mesh'], smooth_node.inputs['Geometry'])
ng.links.new(smooth_node.outputs['Geometry'], set_material_node.inputs['Geometry'])
ng.links.new(set_material_node.outputs['Geometry'], instance_node.inputs['Instance'])
ng.links.new(instance_node.outputs['Instances'], join_geometry_node.inputs['Geometry'])

ng.links.new(sphere_node_1.outputs['Mesh'], smooth_node_1.inputs['Geometry'])
ng.links.new(smooth_node_1.outputs['Geometry'], set_material_node_1.inputs['Geometry'])
ng.links.new(set_material_node_1.outputs['Geometry'], instance_node_1.inputs['Instance'])
ng.links.new(instance_node_1.outputs['Instances'], join_geometry_node.inputs['Geometry'])

ng.links.new(sphere_node_2.outputs['Mesh'], smooth_node_2.inputs['Geometry'])
ng.links.new(smooth_node_2.outputs['Geometry'], set_material_node_2.inputs['Geometry'])
ng.links.new(set_material_node_2.outputs['Geometry'], instance_node_2.inputs['Instance'])
ng.links.new(instance_node_2.outputs['Instances'], join_geometry_node.inputs['Geometry'])

ng.links.new(join_geometry_node.outputs['Geometry'],realize_instances_node.inputs['Geometry'])
ng.links.new(realize_instances_node.outputs['Geometry'],output.inputs['Geometry'])
        #links: ID comparison
ng.links.new(input.outputs['TYPE'], ID_compare_node_1.inputs[0])
ng.links.new(ID_compare_node_1.outputs['Result'], separate_geometry_node_1.inputs['Selection'])
ng.links.new(input.outputs['TYPE'], ID_compare_node_2.inputs[0])
ng.links.new(ID_compare_node_2.outputs['Result'], separate_geometry_node_2.inputs['Selection'])

        #links: Value Color
ng.links.new(instance_node.outputs['Instances'], minmax_node.inputs['Geometry'])
ng.links.new(input.outputs[coloring_field], minmax_node.inputs['Attribute'])
ng.links.new(input.outputs[coloring_field], map_range_node.inputs['Value'])
ng.links.new(minmax_node.outputs['Min'], map_range_node.inputs['From Min'])
ng.links.new(minmax_node.outputs['Max'], map_range_node.inputs['From Max'])
ng.links.new(map_range_node.outputs['Result'], color_node.inputs['Fac'])
ng.links.new(color_node.outputs['Color'], output.inputs['AtomColor'])

    # bonds
print("defining bond geonode modifier")
bpy.context.view_layer.objects.active = bond_obj
bpy.ops.object.modifier_add(type='NODES')
bpy.ops.node.new_geometry_node_group_assign()
ng = bond_obj.modifiers[-1].node_group
nodes = ng.nodes

input = nodes.get('Group Input')
input.location = (-550, 0)
output = nodes.get('Group Output')
output.location = (850, 0)
        # attributes
print("defining bond attributes")
new_attribute = ng.inputs.new('NodeSocketVectorXYZ', "rotation")
bond_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'rotation'
bond_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
new_attribute = ng.inputs.new('NodeSocketFloat', "scale")
bond_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'scale'
bond_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
new_attribute = ng.inputs.new('NodeSocketFloat', "bond_value")
bond_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'bond_value'
bond_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
new_attribute = ng.outputs.new('NodeSocketColor','BondColor')
bond_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'bond_color'
bond_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
        # nodes
print("defining bond nodes")
cyl_node = nodes.new('GeometryNodeMeshCylinder')
cyl_node.location = (-400, 300)
cyl_node.inputs['Radius'].default_value=0.05
cyl_node.inputs['Depth'].default_value=1.0
set_material_node = nodes.new('GeometryNodeSetMaterial')
set_material_node.location = (-50, 300)
set_material_node.inputs['Material'].default_value = bond_value_material
smooth_node = nodes.new('GeometryNodeSetShadeSmooth')
smooth_node.location = (-200, 300)
combine_xyz_node = nodes.new('ShaderNodeCombineXYZ')
combine_xyz_node.location = (-80, 30)
combine_xyz_node.inputs['X'].default_value = 1.0
combine_xyz_node.inputs['Y'].default_value = 1.0
instance_node = nodes.new('GeometryNodeInstanceOnPoints')
instance_node.location = (150, 200)
minmax_node = nodes.new('GeometryNodeAttributeStatistic')
minmax_node.location = (15, -110)
map_range_node = nodes.new('ShaderNodeMapRange')
map_range_node.location = (200, -50)
color_node = nodes.new('ShaderNodeValToRGB')
color_node.location = (400, 0)
realize_instances_node = nodes.new('GeometryNodeRealizeInstances')
realize_instances_node.location = (300, 250)

        #links: geometry
print("defining bond links")
ng.links.new(input.outputs['Geometry'], instance_node.inputs['Points'])
ng.links.new(input.outputs['rotation'], instance_node.inputs['Rotation'])
ng.links.new(input.outputs['scale'], combine_xyz_node.inputs['Z'])
ng.links.new(combine_xyz_node.outputs['Vector'], instance_node.inputs['Scale'])
ng.links.new(cyl_node.outputs['Mesh'], smooth_node.inputs['Geometry'])
ng.links.new(smooth_node.outputs['Geometry'], set_material_node.inputs['Geometry'])
ng.links.new(set_material_node.outputs['Geometry'], instance_node.inputs['Instance'])
ng.links.new(instance_node.outputs['Instances'], realize_instances_node.inputs['Geometry'])
ng.links.new(realize_instances_node.outputs['Geometry'],output.inputs['Geometry'])

        #links: Value Color
ng.links.new(input.outputs['Geometry'], minmax_node.inputs['Geometry'])
ng.links.new(input.outputs['bond_value'], minmax_node.inputs['Attribute'])
ng.links.new(input.outputs['bond_value'], map_range_node.inputs['Value'])
ng.links.new(minmax_node.outputs['Min'], map_range_node.inputs['From Min'])
ng.links.new(minmax_node.outputs['Max'], map_range_node.inputs['From Max'])
ng.links.new(map_range_node.outputs['Result'], color_node.inputs['Fac'])
ng.links.new(color_node.outputs['Color'], output.inputs['BondColor'])

print("defining shader nodes")
# shader nodes
attribute_node = atom_value_material.node_tree.nodes.new("ShaderNodeAttribute")
attribute_node.location = (-280, 300)
attribute_node.attribute_name = "atom_color"
atom_value_material.node_tree.links.new(attribute_node.outputs['Color'], 
                        atom_value_material.node_tree.nodes['Principled BSDF'].inputs['Base Color'])
atom_value_material.node_tree.links.new(attribute_node.outputs['Alpha'], 
                        atom_value_material.node_tree.nodes['Principled BSDF'].inputs['Alpha'])

attribute_node = bond_value_material.node_tree.nodes.new("ShaderNodeAttribute")
attribute_node.attribute_name = "bond_color"
bond_value_material.node_tree.links.new(attribute_node.outputs['Color'], 
                        bond_value_material.node_tree.nodes['Principled BSDF'].inputs['Base Color'])
bond_value_material.node_tree.links.new(attribute_node.outputs['Alpha'], 
                        bond_value_material.node_tree.nodes['Principled BSDF'].inputs['Alpha'])
print("done")

# frame update function
def mesh_update(scene):
    if scene.frame_current <= len(dump_data_list):
        if composite:
            fields, atoms, N, time = dump_data_list[scene.frame_current - 1]
        else:
            fields, atoms, N, time = import_dump.lammps_single(dump_data_list[scene.frame_current - 1])
            bfields, bonds, bN, _ = import_dump.lammps_bond_single(bond_data_list[0], cutoff=2.0, cutoff_index=bond_dist_ix)
        atom_obj = bpy.context.scene.objects['atoms']
        bond_obj = bpy.context.scene.objects['bonds']
        if "ux" in fields:
                x_ix = fields.index("ux")
                y_ix = fields.index("uy")
                z_ix = fields.index("uz")
        else:
                x_ix = fields.index("x")
                y_ix = fields.index("y")
                z_ix = fields.index("z")
        rs = np.array([[atom[x_ix], atom[y_ix], atom[z_ix]] for atom in atoms])
        mesh = bpy.data.meshes.new("mesh")
        mesh.from_pydata(rs,[],[])
        atom_obj.data = mesh
        if np.all(np.arange(1, N + 1) == np.array([a[0] for a in atoms])):
                midpoints = [np.array(atoms[int(np.minimum(b[bond_id1_ix],b[bond_id2_ix])-1)][x_ix:z_ix+1]) - np.array(b[bond_dx_ix:bond_dx_ix+3])/2 for b in bonds]
                rotations = [import_dump.vec2rot(b[bond_dx_ix:bond_dx_ix+3]) for b in bonds]
                scales = [b[bond_dist_ix] for b in bonds]
        else:
                raise Exception("implement bond endpoint finding for non-sequential atom numbering")

        for i, field in enumerate(fields):
                atom_obj.data.attributes.new(name=field.upper(), type='FLOAT', domain='POINT')
                atom_obj.data.attributes[field.upper()].data.foreach_set('value', [atom[i] for atom in atoms])
        mesh = bpy.data.meshes.new("mesh")
        mesh.from_pydata(midpoints,[],[])
        bond_obj.data = mesh

        for i, field in enumerate(bfields):
                bond_obj.data.attributes.new(name=field.upper(), type='FLOAT', domain='POINT')
                bond_obj.data.attributes[field.upper()].data.foreach_set('value', [bond[i] for bond in bonds])
        bond_obj.data.attributes.new(name='rotation', type='FLOAT_VECTOR', domain='POINT')
        rotations = np.array(rotations).reshape(-1)
        bond_obj.data.attributes['rotation'].data.foreach_set('vector', rotations)

        bond_obj.data.attributes.new(name='scale', type='FLOAT', domain='POINT')
        bond_obj.data.attributes['scale'].data.foreach_set('value', scales)


bpy.app.handlers.frame_change_pre.append(mesh_update)
