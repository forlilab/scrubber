from rdkit import Chem
import numpy as np
import math
from .utils import rotation_matrix

def find_n_ring_substituents(mol):
    """
    Find a positively charged N atom in a ring
    with two substituents (non-ring atoms).
    Return (n_idx, sub1_idx, sub2_idx) if found, else None.
    """
    smarts = '[NRX4H1,NRX3]'
    patt = Chem.MolFromSmarts(smarts)

    matches = mol.GetSubstructMatches(patt)

    amines = []

    for match in matches:
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)

        neighbors = n_atom.GetNeighbors()
        ring_neighbors = [a for a in neighbors if a.IsInRing()]
        sub_neighbors = [a for a in neighbors if not a.IsInRing()]

        if len(ring_neighbors) == 2 and len(sub_neighbors) == 2:
            sub1_idx = sub_neighbors[0].GetIdx()
            sub2_idx = sub_neighbors[1].GetIdx()
            amines.append((n_idx, sub1_idx, sub2_idx))
        elif len(ring_neighbors) == 2 and len(sub_neighbors) == 1: # case with 1 substituent
            sub1_idx = sub_neighbors[0].GetIdx()
            amines.append((n_idx, sub1_idx))

    return amines

def get_substituent_tree(mol, root_idx, exclude_idx):
    """
    Return all atom indices in the substituent connected to root_idx,
    excluding the path back to exclude_idx (typically the nitrogen).
    """
    visited = set()
    stack = [root_idx]
    substituent = []

    while stack:
        atom_idx = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        substituent.append(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx != exclude_idx and nbr_idx not in visited:
                stack.append(nbr_idx)

    return substituent

def swap_substituents(mol, coords, n_idx, sub1_idx, sub2_idx, ring_atom_indices):

    """
    Rotate the two substituents around the axis bisecting the angle between the two ring atoms bonded to nitrogen.

    coords: numpy array of shape (N_atoms, 3)
    n_idx: index of nitrogen atom
    sub1_idx: index of first substituent atom
    sub2_idx: index of second substituent atom; may be None
    ring_atom_indices: list of indices of atoms forming the ring

    Returns a copy of the modified coords.
    """
    coords_new = coords.copy()

    # Find the two ring atoms connected to nitrogen
    ring_neighbors = []
    nitrogen = mol.GetAtomWithIdx(n_idx)
    neighbors = nitrogen.GetNeighbors()
    for neighbor in neighbors:
        idx = neighbor.GetIdx()
        if idx not in [sub1_idx, sub2_idx]:
            ring_neighbors.append(neighbor.GetIdx())

    if len(ring_neighbors) != 2:
        raise ValueError("Could not identify two ring atoms bonded to nitrogen")

    r1 = coords[ring_neighbors[0]]
    r2 = coords[ring_neighbors[1]]
    n = coords[n_idx]

    # Define vectors
    v1 = r1 - n
    v2 = r2 - n

    # Axis is bisector of v1 and v2
    bisector = v1 / np.linalg.norm(v1) + v2 / np.linalg.norm(v2)
    bisector /= np.linalg.norm(bisector)

    # Get full substituent trees
    group1 = get_substituent_tree(mol, sub1_idx, n_idx)
    if sub2_idx is not None:
        group2 = get_substituent_tree(mol, sub2_idx, n_idx)
    else:
        group2 = []

    for sub_group in [group1, group2]:
        for atom_idx in sub_group:
            vec = coords[atom_idx] - n
            rotated_vec = np.dot(rotation_matrix(bisector, math.pi), vec)
            coords_new[atom_idx] = n + rotated_vec

    return coords_new


def get_ring_atoms(mol, n_idx):
    """
    Given the nitrogen index, find the ring it's part of and return the ring atom indices.
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if n_idx in ring and (len(ring) == 5 or len(ring) == 6):
            return ring
    return None

# for debbuging .
def molToXYZ(mol, coords):
    """
    Convert an RDKit mol and numpy array of coordinates to XYZ format string.
    """
    lines = [str(mol.GetNumAtoms()), "XYZ generated from RDKit"]
    for atom, pos in zip(mol.GetAtoms(), coords):
        sym = atom.GetSymbol()
        x, y, z = pos
        lines.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines)

# testing
# from rdkit.Chem import SDMolSupplier
# sdf_path = "amine.sdf"  
# suppl = SDMolSupplier(sdf_path, removeHs=False)

# for mol in suppl:
#     if mol is None:
#         continue

#     match = find_n_ring_substituents(mol)
#     if match:
#         n_idx, sub1_idx, sub2_idx = match
#         ring_atoms = get_ring_atoms(mol, n_idx)
#         conf = mol.GetConformer()
#         coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

#         new_coords = swap_substituents(coords, mol, n_idx, sub1_idx, sub2_idx, ring_atoms)

#         print(f"Original coords for substituents: {coords[sub1_idx]}, {coords[sub2_idx]}")
#         print(f"Swapped coords for substituents: {new_coords[sub1_idx]}, {new_coords[sub2_idx]}")
#         break
#     else:
#         print("No matching nitrogen found.")
