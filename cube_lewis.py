import bpy
import bmesh
import numpy as np
from os import scandir
import pyopenvdb as openvdb

def read_cube(path, make_grid=True):
    ldata = np.loadtxt(path)
    xs = ldata[:, 0]
    ys = ldata[:, 1]
    zs = ldata[:, 2]
    rs = ldata[:, 0:3]

    deltas = np.diff(np.sort(xs))
    deltas = deltas[deltas > 0]
    delta = np.min(deltas)

    xs = np.arange(np.min(xs), np.max(xs), delta)
    ys = np.arange(np.min(ys), np.max(ys), delta)
    zs = np.arange(np.min(zs), np.max(zs), delta)

    energies = np.zeros((len(xs), len(ys), len(zs)))
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            for k, z in enumerate(zs):
                eps = 0.001
                mask = (abs(ldata[:, 0] - x) < eps) & (abs(ldata[:, 1] - y) < eps) & (abs(ldata[:, 2] - z) < eps)
                if np.any(mask):
                    energies[i, j, k] = ldata[mask, 3]
                else:
                    energies[i, j, k] = np.nan

    energies = np.nan_to_num(energies, nan=np.max(energies[~np.isnan(energies)]))
    emax = -0.2 #np.max(energies)
    emin = np.min(energies)
    densities = (emax - energies) / (emax - emin)
    densities = np.clip(densities, 0, 1)
    densities = np.nan_to_num(densities, nan=0.0)

    grid = openvdb.FloatGrid()
    grid.copyFromArray(densities.astype(float))
    grid.transform = openvdb.createLinearTransform([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])
    grid.gridClass = openvdb.GridClass.FOG_VOLUME
    grid.name='density'

    return rs, grid


# define input ============================================================================================
cube_location = r"/Users/jd6157/Documents/blender_projects/lewis_densities/coords_energies_2.txt"  # folder if individual. file if composite
atoms_location = r"/Users/jd6157/Documents/blender_projects/lewis_densities/atoms.txt"
# define input ============================================================================================

print("loading cube file")
rs, grid = read_cube(cube_location)
atom_rs = np.loadtxt(atoms_location)

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

# create mesh

#print("creating mesh")
#mesh = bpy.data.meshes.new("mesh")
#points_obj = bpy.data.objects.new('points', mesh)
#bpy.context.scene.collection.children[0].objects.link(points_obj)
#mesh.from_pydata(rs, [], [])
print("creating mesh")
mesh = bpy.data.meshes.new("mesh")
atom_obj = bpy.data.objects.new('atoms', mesh)
bpy.context.scene.collection.children[0].objects.link(atom_obj)
mesh.from_pydata(atom_rs, [], [])
