import bpy
import bmesh
import numpy as np
from os import scandir
import sys
sys.path.insert(0, r"C:\Users\jacks\Desktop\Blender\python")
import import_dump
# define input ============================================================================================
dump_location = r"." # folder if individual. file if composite
composite = False
bond_id1_ix = 0
bond_dx_ix = 2
bond_dist_ix = 5
# define input ============================================================================================

print("loading dump file")
dump_data_list = [f.path for f in scandir(dump_location) if "dump" in f.name and "atom" in f.name]
bond_data_list = [f.path for f in scandir(dump_location) if "dump" in f.name and "bond" in f.name]
bpy.context.scene.frame_end = len(dump_data_list)
dump_data_list.sort(key=import_dump.dumpnum)
bond_data_list.sort(key=import_dump.dumpnum)

# clear current junk handlers
while len(bpy.app.handlers.frame_change_pre) > 0:
    bpy.app.handlers.frame_change_pre.pop()
def mesh_update(scene):
    if 0 < scene.frame_current <= len(dump_data_list):
        if composite:
            fields, atoms, N, time = dump_data_list[scene.frame_current - 1]
        else:
            fields, atoms, N, time = import_dump.lammps_single(dump_data_list[scene.frame_current - 1])
            bfields, bonds, bN, _ = import_dump.lammps_bond_single(bond_data_list[scene.frame_current - 1], cutoff=2.0, cutoff_index=bond_dist_ix)
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
                midpoints = [np.array(atoms[int(np.minimum(b[bond_id1_ix],b[bond_id1_ix+1])-1)][x_ix:z_ix+1]) - np.array(b[bond_dx_ix:bond_dx_ix+3])/2 for b in bonds]
        else:
                atom_id_dict = {i: a[0] for i, a in enumerate(atoms)}
                midpoints = [np.array(atom_id_dict[int(np.minimum(b[bond_id1_ix],b[bond_id1_ix+1]))][x_ix:z_ix+1]) - np.array(b[bond_dx_ix:bond_dx_ix+3])/2 for b in bonds]
        rotations = [import_dump.vec2rot(b[bond_dx_ix:bond_dx_ix+3]) for b in bonds]
        scales = [b[bond_dist_ix] for b in bonds]

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
