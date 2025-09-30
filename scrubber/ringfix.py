import numpy as np
import math
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from typing import List

from rdkit.Chem import AllChem
import random


from .utils import optimize_conformers, add_conformers_to_mol, write_conformers_to_sdf
from .utils import rotation_matrix
from . import amine_flip as am


def norm(v):
    return v / np.sqrt(np.dot(v, v))

def dihedral(A,B,C,D):
    """Calculate dihedral considering A in the beggining"""
    A, B, C, D = [np.array(x) for x in (A,B,C,D)]
    b1 = B - A
    b2 = C - B
    b3 = D - C
    temp = np.linalg.norm(b2) * b1
    y = np.dot(temp, np.cross(b2, b3))
    x = np.dot(np.cross(b1, b2), np.cross(b2, b3))
    return np.arctan2(y, x)

def bite_own_tail_recursively(mol, seed_atom_idx, visited, tail_indices):
    """ find if a substituent leads back to a ring (e.g. in two fused rings)
        returns False for spiro as long as the starting ring atom
        is not included in tail_indices
        populates visited with the atoms that are downstream from a substutient
    """
    does_bite = False
    for neigh in mol.GetAtomWithIdx(seed_atom_idx).GetNeighbors():
        idx = neigh.GetIdx()
        if idx in visited:
            continue
        if idx in tail_indices:
            does_bite = True
            continue
        visited.append(idx)
        does_bite |= bite_own_tail_recursively(mol, idx, visited, tail_indices)
    return does_bite

class RingInfo:
    def __init__(self, coords, indices, debug=False):
        self.upness = {} # up (positive) or down (negative) relative to normal
        self.dot_edges = {} # how well aligned are the edges that flip each atom
        self.normal = {}
        self.centroid = {}
        self.indices = [i for i in indices]

        n = len(indices)
        for i in range(n):
            a = np.array(coords[indices[(i-2) % n]]) # rod1
            b = np.array(coords[indices[(i-1) % n]]) # hinge1
            c = np.array(coords[indices[(i+0) % n]]) # the current corner
            d = np.array(coords[indices[(i+1) % n]]) # hinge2
            e = np.array(coords[indices[(i+2) % n]]) # rod2
            centroid = np.mean([a, b, d, e], axis=0)
            edge1 = norm(b-a)
            edge2 = norm(d-e)
            dot = np.dot(edge1, edge2)

            self.dot_edges[indices[i]] = dot
            self.centroid[indices[i]] = centroid

            # we have four atoms that define two edges, and a centroid.
            # let's iterate over the four possible trianges defined by the centroid
            # and two adjacent atoms, compute the normal for each of the triangles,
            # and check whether the corners are up or down (positive or negative c1/c2)
            upness = 0.
            normals = []
            edge_points = [a, b, d, e]
            for j in range(4):
                v1 = edge_points[(j+1) % 4] - edge_points[j % 4]
                v2 = centroid - edge_points[(j+1) % 4]
                if np.dot(v1, v1) < 1e-3:
                    continue # in 4-member rings, a and e are the same point, skip that triangle
                normal = norm(np.cross(v1, v2))
                upness += np.dot(normal, norm(c - centroid))
                normals.append(normal)
            self.upness[indices[i]] = upness
            self.normal[indices[i]] = norm(np.mean(normals, axis=0))
            if debug:
                print(
                    "corner: %2d," % (indices[i] + 1),
                    "a->b: %2d -> %2d," % (indices[(i-2)%n]+1, indices[(i-1)%n]+1),
                    "e->d: %2d -> %2d," % (indices[(i+2)%n]+1, indices[(i+1)%n]+1),
                    #"angle: %.2f" % np.degrees(np.arccos(dot)),
                    "dot: %.3f" % dot,
                    "upness: %.3f" % upness,
                )


def get_substituents(mol, indices):
    substituents = {}
    for i in range(len(indices)):
        idx = indices[i]
        other_ring_idxs = [j for j in indices if j != idx] # all ring atoms except the current one
        atom_info = {}
        for neigh in mol.GetAtomWithIdx(idx).GetNeighbors():
            info = {
                "nr_neighbors": 0,
                "downstream_indices": set([idx]),
                "bites_own_tail": False,
                "atomic_nr": neigh.GetAtomicNum(),
            }
            neigh_idx = neigh.GetIdx()
            if neigh_idx in other_ring_idxs:
                continue
            if neigh.GetAtomicNum() == 1:
                info["downstream_indices"].add(neigh_idx)
                atom_info[neigh_idx] = info
                continue
            info["nr_neighbors"] = len([a for a in neigh.GetNeighbors() if a.GetAtomicNum() > 1]) - 1
            visited = [idx, neigh_idx]
            info["bites_own_tail"] |= bite_own_tail_recursively(mol, neigh_idx, visited, other_ring_idxs)
            visited.pop(0) # remove ring atom (idx)
            info["downstream_indices"] = info["downstream_indices"].union(visited)
            atom_info[neigh_idx] = info
        substituents[indices[i]] = atom_info
    return substituents

def calc_boat_likeliness(ring_info):
    if len(ring_info.indices) != 6:
        return None # this is only for 6-member rings
    boat_likeliness = 0. # positive: boat-like; negative: chairs
    for i in range(6):
        idx = ring_info.indices[i]
        opposite_corner_idx = ring_info.indices[(i + 3) % 6]
        dot = ring_info.dot_edges[idx]
        u1 = ring_info.upness[idx]
        u2 = ring_info.upness[opposite_corner_idx]
        boat_likeliness += dot * u1 * u2
    return boat_likeliness

def calc_axial_likeliness_old(ringinfo, substituents, coords):
    axial_likeliness = 0.
    for idx in substituents:
        for neigh_idx in substituents[idx]:
            substituent = substituents[idx][neigh_idx]
            if substituent["atomic_nr"] == 1:
                continue
            if substituent["bites_own_tail"]:
                continue
            weight = 1 + substituent["nr_neighbors"]**2 / 3 # Me:1.0, Et:1.3, Ph:2.3, tert-butyl:4.0
            v = norm(coords[neigh_idx] - coords[idx])
            dot = np.dot(v, ringinfo.normal[idx])
            axial_likeliness += weight * abs(dot)
    return axial_likeliness

def calc_anomeric_penalty(mol, substituents, coords):
    pattern = Chem.MolFromSmarts("[*]-[OX2,SX2]-[CX4]-[OX2,SX2,#9,#17,#35,#53]")
    matches = mol.GetSubstructMatches(pattern)
    penalty = 0.
    for i,j,k,l in matches:
        if (
            i in substituents and
            j in substituents and
            k in substituents
        ):
            angle = dihedral(coords[i], coords[j], coords[k], coords[l])
            k1 = -0.5
            k2 = 2.0
            k3 = 1.5
            penalty += k1 * (1.0 + math.cos(1 * angle)) + 1  # guarantee penalty is non-negative
            penalty += k2 * (1.0 + math.cos(2 * angle))
            penalty += k3 * (1.0 + math.cos(3 * angle))
    return penalty

def calc_axial_likeliness(substituents, coords):
    axial_likeliness = 0.
    centroid = np.mean([coords[j] for j in substituents], axis=0)
    normals = []
    ring_idxs = list(substituents.keys())
    for i in range(len(substituents)):
        j = (i + 1) % len(substituents)
        vi = norm(coords[ring_idxs[i]] - centroid)
        vj = norm(coords[ring_idxs[j]] - centroid)
        normals.append(np.cross(vi, vj))
    normal = np.mean(normals, axis=0) # not normalized so length correlates with planarity of ring
    for idx in substituents:
        for neigh_idx in substituents[idx]:
            substituent = substituents[idx][neigh_idx]
            if substituent["atomic_nr"] == 1:
                continue
            weight = 1 + substituent["nr_neighbors"]**2 / 3 # Me:1.0, Et:1.3, Ph:2.3, tert-butyl:4.0
            v = norm(coords[neigh_idx] - coords[idx])
            dot = np.dot(v, normal)
            axial_likeliness += weight * abs(dot)
    return axial_likeliness


def fix_rings(mol: Mol, 
              coords: list, 
              use_energy: bool, 
              energy_threshold: float, 
              debug=False, 
              max_ff_iter: int = 400,
              ff:str = "mmff94"):
    #one_ring_atom_smarts = "[$([R1]),$([R2;x4]);!$([#6;R2;x3]);!$([#6;R1;X3](@=*));!$([#6](=*)(@N))]"
    #smarts = "{s}1{s}{s}{s}{s}{s}1".format(s=one_ring_atom_smarts)

    ring6_smarts = "[*]1[*][*][*][*][*]1"
    amide_smarts = "[NX3]-[CX3]=[O,N,SX1]"
    amide_idxs = mol.GetSubstructMatches(Chem.MolFromSmarts(amide_smarts))
    ring6_rot6_idxs = []
    ring6_rot5_idxs = []
    for idxs in mol.GetSubstructMatches(Chem.MolFromSmarts(ring6_smarts)):
        is_bond_rotatable = []
        for i in range(len(idxs)):
            a = idxs[i]
            b = idxs[(i+1) % len(idxs)]
            bond_type = mol.GetBondBetweenAtoms(a, b).GetBondType()
            is_single_bond = bond_type == Chem.rdchem.BondType.SINGLE
            is_amide_bond = False
            for indices in amide_idxs:
                if (indices[0], indices[1]) == (a, b) or (indices[0], indices[1]) == (b, a):
                    is_amide_bond = True
                    break
            is_bond_rotatable.append(is_single_bond and not is_amide_bond)

        # saturated rings (e.g. cyclohexane, piperidine)
        if sum(is_bond_rotatable) == 6:
            ring6_rot6_idxs.append(idxs)
            coords = convert_boat_to_chair(mol, coords, idxs, debug)

        # 6-member rings with one rigid bond
        if sum(is_bond_rotatable) == 5:
            offset = is_bond_rotatable.index(False) + 1
            # rigid bond between idxs[0] and idxs[-1]
            idxs = [idxs[(offset + i) % 6] for i in range(6)]
            ring6_rot5_idxs.append(idxs)

    # now we run methods that may return more than a single conformation
    coords_list = [coords]
    for idxs in ring6_rot6_idxs:
        tmp = []
        substituents = get_substituents(mol, idxs)
        for coords in coords_list:
            ringinfo = RingInfo(coords, idxs, debug)
            new_coords = expand_reasonable_chairs(coords, idxs, ringinfo, substituents, mol, use_energy, energy_threshold, debug, 0.1, max_ff_iter, ff)
            tmp.extend(new_coords)
        coords_list = tmp
    for idxs in ring6_rot5_idxs:
        tmp = []
        substituents = get_substituents(mol, idxs)
        for coords in coords_list:
            new_coords = expand_ring6_rot5(coords, idxs, substituents)
            tmp.extend(new_coords)
        coords_list = tmp

    return coords_list


def rotate_amine_substituents(mol: Mol, 
                              coords_list: List, 
                              max_ff_iter: int, 
                              energy_threshold: float, 
                              debug: bool, 
                              ff: str = "mmff94"):
    
    """
    Identify amine group in ring structure and swap equatorial and axial
    substituents. Determine optimal structure with energy minimization. 
    """
    
    amine_match = am.find_n_ring_substituents(mol)

    if len(amine_match) == 0:
        return coords_list
    
    for match in amine_match:
        if len(match) == 3:
            n_idx, sub1_idx, sub2_idx = match
        elif len(match) == 2:
            n_idx, sub1_idx = match
            sub2_idx = None
        ring_atoms = am.get_ring_atoms(mol, n_idx)
        tmp = []
        for coords in coords_list:
            
            new_coords =  am.swap_substituents(mol, coords, n_idx, sub1_idx, sub2_idx, ring_atoms)

            mol_with_confs = add_conformers_to_mol(mol, [coords, new_coords])
            # compare new coords with old coords
            optimized_energies = optimize_conformers(mol_with_confs, ff, max_ff_iter)
        
            old_energy = optimized_energies[0][1]
            new_energy = optimized_energies[1][1]

            if debug:
                print(f"Initial energy: {old_energy}")
                print(f"Initial coords:\n {am.molToXYZ(mol, coords)}")
                print(f"Rotated energy: {new_energy}")
                print(f"Rotated coords:\n {am.molToXYZ(mol, new_coords)}")

            # determine optimal structure within energy threshold. 
            if new_energy - old_energy < -energy_threshold:
                tmp.append(new_coords)
            elif new_energy - old_energy > energy_threshold:
                tmp.append(coords)
            else:
                tmp.append(coords)
                tmp.append(new_coords)

        coords_list = tmp

    return coords_list



def expand_ring6_rot5(coords, idxs, substituents, axial_range=0.1, debug=False):
    """
        (rod 1) r1 -- h1 (hinge 1)
                       \
                        c1 (corner 1)
                        |
                        c2 (corner 2)
                       /
       (rod 2) r2 -- h2 (hinge 2)
    """

    input_coords = coords.copy()
    r1, h1, c1, c2, h2, r2 = (coords[idxs[i]].copy() for i in range(6))

    # Consider the h1-c1-c2-h2 four-atom segment. To invert it's dihedral angle,
    # e.g. from +80 to -80, we would rotate around the c1-c2 bond and calculate
    # the new position of h2 (`new_h2` below). Then, we update the position
    # of c1 and c2 by rotating around np.cross(v1, v2) so that the new_h2 point
    # goes back to the original h2 position. Thus, h2 position remains unchanged.
    angle_corners = dihedral(h1, c1, c2, h2)
    new_h2 = np.dot(rotation_matrix(c2-c1, -2*angle_corners), h2-c2) + c2
    v1 = norm(new_h2-h1)
    v2 = norm(h2-h1)
    axis = np.cross(v1, v2)
    angle = np.arccos(np.dot(v1, v2))
    coords -= h1
    coords = rotate_ring_atom(idxs[2], coords, axis, angle, substituents)
    coords = rotate_ring_atom(idxs[3], coords, axis, angle, substituents)

    # Rotate c1 and c2 around h1-h2 axis to restore the middle point
    # between c1 and c2 to the original position
    orig_mid_corner = (c1 + c2) / 2
    new_mid_corner = (coords[idxs[2]] + coords[idxs[3]]) / 2 + h1
    angle = dihedral(orig_mid_corner, h1, h2, new_mid_corner)
    coords = rotate_ring_atom(idxs[2], coords, h1-h2, angle, substituents)
    coords = rotate_ring_atom(idxs[3], coords, h1-h2, angle, substituents)

    # rotate h1 substituents over r1-h1 axis
    angle = dihedral(r2, r1, h1, coords[idxs[2]]) - dihedral(r2, r1, h1, c1)
    coords = rotate_ring_atom(idxs[1], coords, h1-r1, angle, substituents)
    coords += h1

    # rotate h2 substituents over r2-h2 axis
    angle = dihedral(r1, r2, h2, coords[idxs[3]]) - dihedral(r1, r2, h2, c2)
    coords = rotate_ring_atom(idxs[4], coords-h2, h2-r2, angle, substituents) + h2

    # return coords with least axial likeliness
    input_axial = calc_axial_likeliness(substituents, input_coords)
    mod_axial = calc_axial_likeliness(substituents, coords)
    delta_axial = mod_axial - input_axial
    if delta_axial < -axial_range:
        return [coords]
    elif delta_axial > axial_range:
        return [input_coords]
    else:
        return [input_coords, coords]


def expand_reasonable_chairs(
        coords: list, 
        idxs: tuple, 
        ringinfo: RingInfo, 
        substituents: dict, 
        mol: Mol, 
        use_energy: bool, 
        energy_threshold: float, 
        debug: bool, 
        axial_likeliness_range=0.1,
        max_ff_iter = 400, 
        ff: str = "mmff94"):
    

    # check for amine flip in initially supplied coords
    c = rotate_amine_substituents(mol, [coords], max_ff_iter, energy_threshold, debug, ff)
    coords = c[0]

    if len(idxs) != 6:
        raise RuntimeError("length of idxs is %d but must be 6" % (len(idxs)))
    if calc_boat_likeliness(ringinfo) >= -2:
        return [coords]

    # all corners rotate in a chair-to-other-chair flip, bail out if any
    # substituent leads back to a different ring atom ("bites own tail")
    for idx in substituents:
        for neigh_idx in substituents[idx]:
            if substituents[idx][neigh_idx]["bites_own_tail"]:
                return [coords]

    starting_axial_likeliness = calc_axial_likeliness(substituents, coords)
    starting_axial_likeliness += calc_anomeric_penalty(mol, substituents, coords)

    # find best corner to flip based on alignment of rod-hinge vectors
    # as well as the number of substituents on both corners
    highest_flip_score = -1e10
    for i in range(3):
        c1 = idxs[i]
        c2 = idxs[(i+3) % 6]
        dot = ringinfo.dot_edges[idxs[i]] # co-linearity of rod1-hinge1 and rod2-hinge2 axis
        substituent_weight = 0.
        for idx in [c1, c2]:
            for neigh_idx in substituents[idx]:
                sw = 1 + substituents[idx][neigh_idx]["nr_neighbors"]**2 / 3
                substituent_weight -= 0.05 * sw # low priority relative to rod-hinge alignment
        score = dot**2 + substituent_weight
        if score > highest_flip_score:
            best_index = i
            highest_flip_score = score

    a = coords[idxs[(best_index - 2) % 6]]
    b = coords[idxs[(best_index - 1) % 6]]
    c = coords[idxs[(best_index + 0) % 6]]
    d = coords[idxs[(best_index + 1) % 6]]
    e = coords[idxs[(best_index + 2) % 6]]
    c2= coords[idxs[(best_index + 3) % 6]]
    centroid = ringinfo.centroid[idxs[best_index]]
    rotangle  = 2. * (math.pi-dihedral(centroid, b, d, c))
    rotangle2 = 2. * (math.pi-dihedral(centroid, e, a, c2))
    newpos = rotate_corner(idxs[best_index],           ringinfo, substituents, coords, rotangle)
    newpos = rotate_corner(idxs[(best_index + 3) % 6], ringinfo, substituents, newpos, rotangle2)
    new_axial_likeliness = calc_axial_likeliness(substituents, newpos)
    new_axial_likeliness += calc_anomeric_penalty(mol, substituents, newpos)

    # now check for amine flip on new coordinates
    c = rotate_amine_substituents(mol, [newpos], max_ff_iter, energy_threshold, debug, ff)
    newpos = c[0]

    # a little more debugging
    if debug:
        print("Mol 1")
        print(am.molToXYZ(mol, coords))
        print("Mol 2")
        print(am.molToXYZ(mol, newpos))


    if (use_energy):
        ## calculate correct conformation by energy comparison ######
        ## MMFF94 forcefield
        
        if debug:
            print("Optimizing ring geometries")

        mol_with_confs = add_conformers_to_mol(mol, [coords, newpos])
        if debug: 
            oldmol = add_conformers_to_mol(mol, [coords])
            newmol = add_conformers_to_mol(mol, [newpos])

        optimized_energies = optimize_conformers(mol_with_confs, ff, max_ff_iter)

        # optimized_energy_old = optimize_conformers(oldmol, ff)
        # optimized_energy_new = optimize_conformers(newmol, ff)

        # Print the optimized energy values
        old_energy = optimized_energies[0][1]
        new_energy = optimized_energies[1][1]
        # old_energy = optimized_energy_old[0][1]
        # new_energy = optimized_energy_new[0][1]

        if debug:
            print(f"Old Energy = {old_energy:.4f} kcal/mol")
            print(f"New Energy = {new_energy:.4f} kcal/mol")
            # print(f"Number of conformers: {mol_with_confs.GetNumConformers()}")
            write_conformers_to_sdf(oldmol, f"rings-1-test.sdf")
            write_conformers_to_sdf(newmol, f"rings-2-test.sdf")
            # for conf_id, energy in optimized_energies:
            #     print(f"Conformer {conf_id}: Energy = {energy:.4f} kcal/mol")
                
        if new_energy - old_energy < -energy_threshold:
            return [newpos]
        elif new_energy - old_energy > energy_threshold:
            return [coords]
        else:
            return [coords, newpos]

    #####################################################
    else: 
        delta_axial_likeliness = new_axial_likeliness - starting_axial_likeliness
        if starting_axial_likeliness < 0.001 and new_axial_likeliness < 0.001: # no subs?
            return [coords] # avoids expanding nr confs when unnecessary
        elif delta_axial_likeliness > axial_likeliness_range:
            return [coords]
        elif delta_axial_likeliness < -axial_likeliness_range:
            return [newpos]
        else:
            return [coords, newpos] # new and starting similar, return both


def convert_boat_to_chair(mol, coords, idxs, debug):


    ringinfo = RingInfo(coords, idxs, debug)
    if calc_boat_likeliness(ringinfo) < -2: # chair already
        return coords

    substituents = get_substituents(mol, idxs)
    scores = []
    for i, idx in enumerate(idxs):
        opposite_corner_idx = idxs[(i + 3) % 6]
        a_idx = idxs[(i - 2) % 6] # hinge1
        b_idx = idxs[(i - 1) % 6] # hinge1
        d_idx = idxs[(i + 1) % 6] # hinge2
        e_idx = idxs[(i + 2) % 6] # hinge1
        opposite_upness = ringinfo.upness[opposite_corner_idx]
        this_upness = ringinfo.upness[idx]
        # About flip_score:
        # If a corner is nearly in plane, and the other is far from plane,
        # we want to flip the one that is nearly in plane, thus we take
        # the modulo of the magnitude of the other corner.
        # The product of `this_upness` by `opposite_upness` is negative for
        # chairs and positive for boats.
        flip_score = abs(opposite_upness) * this_upness * opposite_upness
        flip_score *= ringinfo.dot_edges[idx]**2
        for k in [idx, b_idx, d_idx]:
            data = substituents[k]
            for subst_idx in data:
                info = data[subst_idx]
                if debug:
                    print("corner: %2d, subst: %2d, nr_neigh: %2d, bites own tail: %s" % (k+1, subst_idx+1, info["nr_neighbors"], info["bites_own_tail"]))
                weight = 1 + info["nr_neighbors"]**2 / 3 # Me:1.0, Et:1.3, Ph:2.3, tert-butyl:4.0
                flip_score -= 0.0 * weight + 1000 * info["bites_own_tail"]
        scores.append(flip_score)
        if debug:
            print("corner: %d, flip_score: %.2f, this_up=%.2f, opposite_up=%.2f, dot_edges=%.2f" % (idx+1, flip_score, this_upness, opposite_upness, ringinfo.dot_edges[idx]))

    if debug:
        print("Starting boat_likeliness = %.3f" % calc_boat_likeliness(ringinfo))
        print("Starting axial_score = %.3f" % calc_axial_likeliness(substituents, coords))

    best_coords = coords.copy()
    best_score = float("inf")
    for i, idx in enumerate(idxs):
        flip_score = scores[i]
        if flip_score < 0: # some mildly negative flip_scores may benefit from rotation
            continue

        angle = dihedral(ringinfo.centroid[idx], coords[b_idx], coords[d_idx], coords[idx])
        angle2 = dihedral(ringinfo.centroid[idx], coords[e_idx], coords[a_idx], coords[opposite_corner_idx])
        if debug:
            if angle * angle2 < 0:
                print('WARNING: expected angles of opposing corners to be of same signs')

        # if the starting angle is large (less than 150 deg) rotate to the
        # inverse of that. Otherwise rotate to 150 deg (5 * pi / 6)
        target_angle=3*math.pi/4
        target_magnitude = min(abs(angle), target_angle)
        rotangle = math.pi - abs(angle) + math.pi - target_magnitude
        if angle2 < 0:
            rotangle *= -1

        newpos = rotate_corner(idx, ringinfo, substituents, coords, rotangle)
        new_info = RingInfo(newpos, idxs)
        new_boat_likeliness = calc_boat_likeliness(new_info)
        new_axial_likeliness = calc_axial_likeliness(substituents, newpos)
        new_anomeric_penalty = calc_anomeric_penalty(mol, substituents, newpos)
        score = new_boat_likeliness + new_axial_likeliness + new_anomeric_penalty
        if score < best_score:
            best_coords = newpos
            best_score = score
        if debug:
            print("Boatifying atom %2d" % (idx+1), flip_score)
            print("new boat_likeliness = %.3f" % new_boat_likeliness)
            print("new axial_score = %.3f" % new_axial_likeliness)

    return best_coords

def rotate_corner(idx, ringinfo, substituents, coords, rotangle):
    coords = coords.copy()
    idx = ringinfo.indices.index(idx)
    c_idx = ringinfo.indices[idx]            # corner
    b_idx = ringinfo.indices[(idx - 1) % 6]  # hinge1
    d_idx = ringinfo.indices[(idx + 1) % 6]  # hinge2
    c2_idx = ringinfo.indices[(idx + 3) % 6]  # hinge2
    c = coords[c_idx].copy() # this corner
    a = coords[ringinfo.indices[(idx - 2) % 6]].copy() # rod1
    e = coords[ringinfo.indices[(idx + 2) % 6]].copy() # rod2
    b = coords[b_idx].copy() # hinge1
    d = coords[d_idx].copy() # hinge2
    c2 = coords[c2_idx].copy() # opposite corner

    coords = rotate_ring_atom(c_idx, coords-b, d - b, rotangle, substituents)
    newcorner = coords[c_idx] + b

    # rotate hinge1
    rotangle = dihedral(c2, a, b, newcorner) - dihedral(c2, a, b, c)
    coords = rotate_ring_atom(b_idx, coords, b - a, rotangle, substituents)

    # rotate hinge2
    rotangle = dihedral(c2, e, d, newcorner) - dihedral(c2, e, d, c)
    coords = rotate_ring_atom(d_idx, coords-d+b, d - e, rotangle, substituents)
    coords += d
    return coords

def rotate_ring_atom(index, coords, rotaxis, angle, substituents):
    affected = set([index])
    for _, info in substituents[index].items():
        affected = affected.union(info["downstream_indices"])
    for i in affected:
        coords[i] = np.dot(rotation_matrix(rotaxis, angle), coords[i])
    return coords


