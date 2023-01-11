import numpy as np
import math
from rdkit import Chem

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
        self.angle = {} # rotation state around the axis that rotates each atom
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
            self.angle[indices[i]] = dihedral(centroid, b, d, c)
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

def calc_axial_likeliness(ringinfo, substituents, coords):
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
    
def convert_boats_to_chairs(mol, coords, debug=False):
    one_ring_atom_smarts = "[$([R1]),$([R2;x4]);!$([#6;R2;x3]);!$([#6;R1;X3](@=*));!$([#6](=*)(@N))]"
    smarts = "{s}1{s}{s}{s}{s}{s}1".format(s=one_ring_atom_smarts)
    coords = coords.copy()
    for idxs in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
        coords = convert_boat_to_chair(mol, coords, idxs, debug)
    return coords
        

def convert_boat_to_chair(mol, coords, idxs, debug):

    ringinfo = RingInfo(coords, idxs, debug)
    substituents = get_substituents(mol, idxs)
    scores = []
    for i, idx in enumerate(idxs):
        opposite_corner_idx = idxs[(i + 3) % 6]
        b_idx = idxs[(i - 1) % 6] # hinge1
        d_idx = idxs[(i + 1) % 6] # hinge2
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
        print("Starting axial_score = %.3f" % calc_axial_likeliness(ringinfo, substituents, coords))

    best_coords = coords.copy()
    best_score = float("inf")
    for i, idx in enumerate(idxs):
        flip_score = scores[i]
        if flip_score < 0: # some mildly negative flip_scores may benefit from rotation
            continue
    
        newpos = unboatify_atom(idx, ringinfo, substituents, coords)
        new_info = RingInfo(newpos, idxs)
        new_boat_likeliness = calc_boat_likeliness(new_info)
        new_axial_likeliness = calc_axial_likeliness(new_info, substituents, newpos)
        score = new_boat_likeliness + new_axial_likeliness
        if score < best_score:
            best_coords = newpos
            best_score = score
        if debug:
            print("Boatifying atom %2d" % (idx+1), flip_score)
            print("new boat_likeliness = %.3f" % new_boat_likeliness)
            print("new axial_score = %.3f" % new_axial_likeliness)

    return best_coords

def unboatify_atom(idx, ringinfo, substituents, coords, target_angle=3*math.pi/4):
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

    angle = ringinfo.angle[c_idx]
    angle2 = ringinfo.angle[c2_idx]
    if angle * angle2 < 0:
        print('WARNING: expected angles of opposing corners to be of same signs')
    
    # if the starting angle is large (less than 150 deg) rotate to the
    # inverse of that. Otherwise rotate to 150 deg (5 * pi / 6)
    target_magnitude = min(abs(angle), target_angle)
    rotangle = math.pi - abs(angle) + math.pi - target_magnitude
    if angle2 < 0:
        rotangle *= -1
    hinge_axis = d - b
    coords -= b
    affected = set([c_idx])
    for _, info in substituents[c_idx].items():
        affected = affected.union(info["downstream_indices"])
    for i in affected:
        coords[i] = np.dot(rotation_matrix(hinge_axis, rotangle), coords[i])
    newcorner = coords[c_idx] + b

    # rotate hinge1
    rotangle = dihedral(c2, a, b, newcorner) - dihedral(c2, a, b, c)
    affected = set([b_idx])
    for _, info in substituents[b_idx].items():
        affected = affected.union(info["downstream_indices"])
    for i in affected:
        coords[i] = np.dot(rotation_matrix(b - a, rotangle), coords[i])
    coords += b - d

    # rotate hinge2
    rotangle = dihedral(c2, e, d, newcorner) - dihedral(c2, e, d, c)
    affected = set([d_idx])
    for _, info in substituents[d_idx].items():
        affected = affected.union(info["downstream_indices"])
    for i in affected:
        coords[i] = np.dot(rotation_matrix(d - e, rotangle), coords[i])
    coords += d
    return coords


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    source: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """

    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
