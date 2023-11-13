import bpy
import bmesh
import numpy as np
from os import scandir
import sys
sys.path.insert(0, r"C:\Users\jacks\Desktop\Blender\python")
import import_dump

# define input ============================================================================================
dump_location = r"test_dumps" # folder if individual. file if composite
composite = False
# define input ============================================================================================

print("loading dump file")
dumpfiles = [f.path for f in scandir(dump_location) if "dump" in f.name]
bpy.context.scene.frame_end = len(dumpfiles)
dumpfiles.sort(key=import_dump.dumpnum)
if composite:
    dump_data_list = import_dump.lammps_composite(dump_location)

# clear current junk handlers
while len(bpy.app.handlers.frame_change_pre) > 0:
    bpy.app.handlers.frame_change_pre.pop()


def mesh_update(scene):
    if scene.frame_current <= len(dumpfiles):
        if composite:
            fields, atoms, N, time = dump_data_list[scene.frame_current - 1]              
        else:
            fields, atoms, N, time = import_dump.lammps_single(dumpfiles[scene.frame_current - 1])

        atom_obj = bpy.context.scene.objects['atoms']
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
        for i, field in enumerate(fields):
                atom_obj.data.attributes.new(name=field.upper(), type='FLOAT', domain='POINT')
                atom_obj.data.attributes[field.upper()].data.foreach_set('value', [atom[i] for atom in atoms])


bpy.app.handlers.frame_change_pre.append(mesh_update)