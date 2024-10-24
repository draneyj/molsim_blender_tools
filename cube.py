import bpy
import bmesh
import numpy as np
from os import scandir
import pyopenvdb as openvdb

def read_cube(path, make_grid=True):
    f = open(path, "r")
    lines = f.readlines()
    lines = lines[2:]
    N_atoms = int(lines[0].split()[0])
    N_points_x = int(lines[1].split()[0])
    x_step = float(lines[1].split()[1])
    N_points_y = int(lines[2].split()[0])
    y_step = float(lines[2].split()[2])
    N_points_z = int(lines[3].split()[0])
    z_step = float(lines[3].split()[3])
    atom_lines = lines[4:4+N_atoms]
    rs = np.array([[float(x) for x in line.split()[-3:]] for line in atom_lines])
    ids = np.array([int(line.split()[0]) for line in atom_lines])
    qs = np.array([float(line.split()[1]) for line in atom_lines])
    if make_grid:
        density_lines = lines[4+N_atoms:]
        densities = np.zeros((N_points_x, N_points_y, N_points_z))
        i = j = k = 0
        for line in density_lines:
            for chunk in line.split():
                densities[i, j, k] = float(chunk)
                k += 1
                if k == N_points_z:
                    k = 0
                    j += 1
                    if j == N_points_y:
                        j = 0
                        i += 1
        
        dmax = np.max(densities)
        dmin = np.min(densities)
        densities = (densities - dmin) / (dmax - dmin)
        densities = np.clip(densities, 0, 1)
        
        grid = openvdb.FloatGrid()
        grid.copyFromArray(densities.astype(float))
        grid.transform = openvdb.createLinearTransform([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])
        grid.gridClass = openvdb.GridClass.FOG_VOLUME
        grid.name='density'
    else:
        grid=None
    f.close()
    return rs, ids, qs, grid


# define input ============================================================================================
cube_location = r"/Users/jd6157/Documents/blender_projects/cl2diss/cubes"  # folder if individual. file if composite
composite = False
coloring_field = 'q'.upper()  # name of heading to color by. Can change later in geonodes. must be uppercase
# define input ============================================================================================

print("loading cube file")
cube_data_list = [f.path for f in scandir(cube_location) if "cube" in f.name]
cube_data_list.sort()
fields = ["x", "y", "z", "id", "q"]
rs, ids, qs, grid = read_cube(cube_data_list[0])

bpy.context.scene.frame_end = len(cube_data_list)
# define meshes
print("creating meshes")

mesh = bpy.data.meshes.new("mesh")
atom_obj = bpy.data.objects.new('atoms', mesh)
bpy.context.scene.collection.children[0].objects.link(atom_obj)

# create mesh geometry data
print("defining geometry data")
mesh.from_pydata(rs, [], [])

print("defining coloring data")
# define coloring attribute data
for field_name, field in zip(["id", "q"], [ids, qs]):
        atom_obj.data.attributes.new(name=field_name.upper(), type='FLOAT', domain='POINT')
        atom_obj.data.attributes[field_name.upper()].data.foreach_set('value', field)

print("defining materials")

# make valued materials
atom_value_material = bpy.data.materials.new(name="atom_val_mat")
atom_value_material.use_nodes = True

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
new_attribute = ng.interface.new_socket(socket_type="NodeSocketFloat", name="TYPE", in_out="INPUT")
atom_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'TYPE'
atom_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
coloring_attribute_socket_name = "Coloring Attribute"
new_attribute = ng.interface.new_socket(socket_type='NodeSocketFloat', name=coloring_attribute_socket_name, in_out="INPUT")
atom_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = coloring_field
atom_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
new_attribute = ng.interface.new_socket(socket_type='NodeSocketColor', name='AtomColor', in_out="OUTPUT")
atom_obj.modifiers[-1][new_attribute.identifier+'_attribute_name'] = 'atom_color'
atom_obj.modifiers[-1][new_attribute.identifier+'_use_attribute'] = True
        # nodes
print("defining atom nodes")        
sphere_node = nodes.new('GeometryNodeMeshUVSphere')
sphere_node.location = (-1000, -600)
sphere_node.inputs['Radius'].default_value=0.4
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
ng.links.new(input.outputs[coloring_attribute_socket_name], minmax_node.inputs['Attribute'])
ng.links.new(input.outputs[coloring_attribute_socket_name], map_range_node.inputs['Value'])
ng.links.new(minmax_node.outputs['Min'], map_range_node.inputs['From Min'])
ng.links.new(minmax_node.outputs['Max'], map_range_node.inputs['From Max'])
ng.links.new(map_range_node.outputs['Result'], color_node.inputs['Fac'])
ng.links.new(color_node.outputs['Color'], output.inputs['AtomColor'])

print("defining shader nodes")
# shader nodes
attribute_node = atom_value_material.node_tree.nodes.new("ShaderNodeAttribute")
attribute_node.location = (-280, 300)
attribute_node.attribute_name = "atom_color"
atom_value_material.node_tree.links.new(attribute_node.outputs['Color'], 
                        atom_value_material.node_tree.nodes['Principled BSDF'].inputs['Base Color'])
atom_value_material.node_tree.links.new(attribute_node.outputs['Alpha'], 
                        atom_value_material.node_tree.nodes['Principled BSDF'].inputs['Alpha'])

print("done")

# Volume ============================================================================================
print("creating volume")
vscale = 0.25
bpy.ops.object.volume_add(align='WORLD', location=(0, 0, 0), scale=(vscale, vscale, vscale))
volume_obj = bpy.context.scene.objects['Volume']
openvdb.write("/tmp/volume.vdb", grid)
volume_obj.data.filepath = "/tmp/volume.vdb"
volume_material = bpy.data.materials.new(name="volume_mat")
volume_material.use_nodes = True
volume_obj.data.materials.append(volume_material)

# frame update function
def mesh_volume_update(scene):
    if scene.frame_current <= len(cube_data_list):
        print("loading cube file")
        vdb_path = os.path.join(cube_location,f"{scene.frame_current:04d}.vdb")
        make_cube = not os.path.exists(vdb_path)
        rs, ids, qs, grid = read_cube(cube_data_list[scene.frame_current - 1], make_cube)

        atom_obj = bpy.context.scene.objects['atoms']
        mesh = bpy.data.meshes.new("mesh")
        mesh.from_pydata(rs,[],[])
        atom_obj.data = mesh
        for field_name, field in zip(["id", "q"], [ids, qs]):
                atom_obj.data.attributes.new(name=field_name.upper(), type='FLOAT', domain='POINT')
                atom_obj.data.attributes[field_name.upper()].data.foreach_set('value', field)
        # os.remove("/tmp/volume.vdb")
        if make_cube:
            openvdb.write(vdb_path, grid)
        volume_obj = bpy.context.scene.objects['Volume']
        volume_obj.data.filepath = vdb_path


bpy.app.handlers.frame_change_pre.append(mesh_volume_update)