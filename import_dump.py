import numpy as np
from os.path import basename
from scipy.spatial.transform import Rotation as R

def lammps_composite(path, idmap = {1:6,2:18}, loadingbar=True):
    dump_data_list = []
    dimensions = np.zeros((3,2))
    N = 0
    time = 0.0
    reading_atoms = False
    with open(path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break  # stop at eof

            if "ITEM" in line:
                if "TIMESTEP" in line:
                    line = f.readline()
                    time = int(line)
                elif "NUMBER OF ATOMS" in line:
                    line = f.readline()
                    N = int(line)
                elif "BOX BOUNDS" in line:
                    for j in range(3):
                        line = f.readline()
                        chunks = line.split()
                        dimensions[j,0] = float(chunks[0])
                        dimensions[j,1] = float(chunks[1])
                elif "ATOMS" in line:
                    reading_atoms = True
                    fields = line.split()[2:]
                    print(fields)
            if reading_atoms:
                atoms = [np.zeros(len(fields))] * N
                for atom_i in range(N):
                    line = f.readline()
                    atoms[atom_i] = [float(chunk) for chunk in line.split()]
                reading_atoms = False
                dump_data_list.append((fields, atoms, N, time))
    return dump_data_list

def lammps_single(path):
    dimensions = np.zeros((3,2))
    N = 0
    time = 0.0
    reading_atoms = False
    with open(path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break  # stop at eof

            if "ITEM" in line:
                if "TIMESTEP" in line:
                    line = f.readline()
                    time = int(line)
                elif "NUMBER OF ATOMS" in line:
                    line = f.readline()
                    N = int(line)
                elif "BOX BOUNDS" in line:
                    for j in range(3):
                        line = f.readline()
                        chunks = line.split()
                        dimensions[j,0] = float(chunks[0])
                        dimensions[j,1] = float(chunks[1])
                elif "ATOMS" in line:
                    reading_atoms = True
                    fields = line.split()[2:]
                    print(fields)
            if reading_atoms:
                atoms = [np.zeros(len(fields))] * N
                for atom_i in range(N):
                    line = f.readline()
                    atoms[atom_i] = [float(chunk) for chunk in line.split()]
                reading_atoms = False
    
    return fields, atoms, N, time

def lammps_bond_single(path, cutoff_index=0, cutoff=False):
    dimensions = np.zeros((3,2))
    N = 0
    time = 0.0
    reading_entries = False
    with open(path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break  # stop at eof

            if "ITEM" in line:
                if "TIMESTEP" in line:
                    line = f.readline()
                    time = int(line)
                elif "NUMBER OF ENTRIES" in line:
                    line = f.readline()
                    N = int(line)
                elif "BOX BOUNDS" in line:
                    for j in range(3):
                        line = f.readline()
                        chunks = line.split()
                        dimensions[j,0] = float(chunks[0])
                        dimensions[j,1] = float(chunks[1])
                elif "ENTRIES" in line:
                    reading_entries = True
                    fields = line.split()[2:]
                    print(fields)
            if reading_entries:
                bonds = [np.zeros(len(fields))] * N
                bond_i = 0
                for i in range(N):
                    line = f.readline()
                    bonds[bond_i] = [float(chunk) for chunk in line.split()]
                    if not (cutoff and bonds[bond_i][cutoff_index] > cutoff):
                        bond_i += 1
                bonds = bonds[:bond_i]
                reading_entries = False
    
    return fields, bonds, bond_i, time


def dumpnum(path):
    if "_" in basename(path):
        return dumpnums(path)
    return float(basename(path).replace("dump", "").replace("atoms", "").replace("bonds", "")[:-1])

def dumpnums(path):
    nums = basename(path).replace("dump", "").replace("atoms", "").replace("bonds", "")[:-1].split("_")
    nums = [float(n) for n in nums]
    return nums[0] * 1e8 + nums[1]  # hack for now since I didn't sleep. LOOK HERE IF THERE IS A BUG ON ORDERING

def dumpnum2fields_c(path):
    nums = [n for n in basename(path).replace("dump", "").split(".") if len(n)>0]
    return int(nums[0])

def dumpnum2fields_step(path):
    nums = [n for n in basename(path).replace("dump", "").split(".") if len(n)>0]
    if len(nums) == 1:
        return 0
    return int(nums[1])

def vec2rot(v):
    n = np.linalg.norm(v)
    r = R.from_matrix([[0, 0, v[i]/n] for i in range(3)])
    return r.as_euler('xyz', degrees=False)