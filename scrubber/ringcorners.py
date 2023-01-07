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

def process_ring6(coords, indices, debug=False):
    assert(len(indices) == 6)
    # for all three pairs of opposing edges, compute their angle

    boat_score = 0.0
    flip_score = [None] * len(indices) # higher score means we should flip atom to convert boat to chair
    dot_edges = [None] * len(indices) # how well aligned are the edges that flip each atom
    angles = [None] * len(indices)
    centroids = [None] * len(indices)

    for i in range(3):
        a = np.array(coords[indices[(i+0) % 6]]) # rod1
        b = np.array(coords[indices[(i+1) % 6]]) # hinge1
        c = np.array(coords[indices[(i+3) % 6]]) # hinge2
        d = np.array(coords[indices[(i+4) % 6]]) # rod2
        corner1 = coords[indices[(i+2) % 6]]
        corner2 = coords[indices[(i+5) % 6]]
        centroid = np.mean([a, b, c, d], axis=0)
        edge1 = norm(b-a)
        edge2 = norm(c-d)
        dot = np.dot(edge1, edge2)
        dot_edges[(i+2) % 6] = dot
        dot_edges[(i+5) % 6] = dot
        angles[(i + 2) % 6] = dihedral(centroid, b, c, corner1)
        angles[(i + 5) % 6] = dihedral(centroid, d, a, corner2)
        #centroids[(i + 2) % 6] = centroid
        #centroids[(i + 5) % 6] = centroid

        if debug:
            print(
                "edge1: %2d %2d" % (indices[(i+0)%6]+1, indices[(i+1)%6]+1),
                "edge2: %2d %2d" % (indices[(i+3)%6]+1, indices[(i+4)%6]+1),
                "angle: %.2f" % np.degrees(np.arccos(dot)),
                "dot: %.3f" % dot,
            )

        # up / down
        edge_points = [a, b, c, d]
        # c1 and c2 are the out of plane magnitudes. They can be either
        # positive or negative depending on direction in which the atom
        # projects out of plane.
        c1 = 0.0
        c2 = 0.0
        # we have four atoms that define two edges, and a centroid.
        # let's iterate over the four possible trianges defined by the centroid
        # and two adjacent atoms, compute the normal for each of the triangles,
        # and check whether the corners are up or down (positive or negative c1/c2)
        #### normals = []
        for j in range(4):
            v1 = norm(edge_points[(j+1) % 4] - edge_points[j % 4])
            v2 = norm(centroid - edge_points[(j+1) % 4])
            normal = np.cross(v1, v2)
            c1 += np.dot(normal, norm(corner1-centroid))
            c2 += np.dot(normal, norm(corner2-centroid))
            ####### normals.append(normal)

        # if both corners are up or down, c1 * c2 is positive
        # if one corner is up and the other is down, c1 * c2 is negative
        # we multiply by dot to weight by the quality of the underlying
        # plane
        boat_score += dot * c1 * c2
        # if a corner is nearly in plane, and the other is far from plane,
        # we want to flip the one that is nearly in plane, thus we take
        # the modulo of the magnitude of the other corner
        flip_score[(i+2) % 6] = abs(c2) * c1 * c2
        flip_score[(i+5) % 6] = abs(c1) * c1 * c2

        #####avg_normal = norm(np.mean(normals, axis=0))
        #####orientations[(i+2) % 6] = (avg_normal, c2) 
        #####orientations[(i+5) % 6] = (avg_normal, c1)
        if debug:
            print("c1: %.3f" % c1, "c2: %.3f" % c2, "product = %.3f" % (c1*c2))
    return boat_score, flip_score, dot_edges, angles#, centroids
    #print("boat score: %.2f" % boat_score)

def bite_own_tail_recursively(mol, seed_atom_idx, visited, tail_indices):
    does_bite = False
    for neigh in mol.GetAtomWithIdx(seed_atom_idx).GetNeighbors():
        idx = neigh.GetIdx()
        print("entering with %d, looping over %d" % (seed_atom_idx, idx), visited)
        if idx in visited:
            continue
        if idx in tail_indices:
            does_bite = True
            continue
        print("here")
        visited.append(idx)
        does_bite |= bite_own_tail_recursively(mol, idx, visited, tail_indices)
    return does_bite
    
def get_substituents(mol, indices, debug=False):
    data = {}
    for idx in indices:
        other_ring_idxs = [i for i in indices if i != idx] # all ring atoms except the current one
        atom_info = {}
        for neigh in mol.GetAtomWithIdx(idx).GetNeighbors():
            info = {
                "nr_neighbors": 0,
                "downstream_indices": set([idx]),
                "bites_own_tail": False,
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
        data[idx] = atom_info
    return data

def convert_boats_to_chairs(mol, coords):
    smarts = "[$([R1]),$([R2;x4]);!$([#6;R2;x3]);!$([#6;R1;X3](@=*));!$([#6](=*)(@N))]" # one of six atoms
    smarts = smarts + "1" + smarts + smarts + smarts + smarts + smarts + "1"
    #if coords is None:
    #    coords = mol.GetConformer().GetPositions()
    new_coords = coords.copy()
    for idxs in mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)):
        score, flip_score, dot_edges, angles = process_ring6(coords, idxs, debug=True)
        data = get_substituents(mol, idxs)
        best_score = float("-inf")
        for i, idx in enumerate(idxs):
            score = flip_score[i] * dot_edges[i]**2
            for subst_idx in data[idx]:
                info = data[idx][subst_idx]
                score -= 0.1 + 0.1 * info["nr_neighbors"] + 1000 * info["bites_own_tail"]
            if score > best_score:
                best_score = score
                best_index = i
        if best_score > 0: # there's probably room for improvment over this simplyfied "if > 0"
            corner_idx = idxs[best_index]
            hinge1_idx = idxs[(best_index - 1) % 6]
            hinge2_idx = idxs[(best_index + 1) % 6]
            corner  = new_coords[corner_idx].copy()
            rod1    = new_coords[idxs[(best_index - 2) % 6]].copy()
            rod2    = new_coords[idxs[(best_index + 2) % 6]].copy()
            hinge1  = new_coords[hinge1_idx].copy()
            hinge2  = new_coords[hinge2_idx].copy()
            ocorner = new_coords[idxs[(best_index + 3) % 6]].copy()

            angle = angles[best_index]
            angle2 = angles[(best_index + 3) % 6] # the opposite corner
            if angle * angle2 < 0:
                print('WARNING: expected angles of opposing corners to be of same signs')
            
            # if the starting angle is large (less than 150 deg) rotate to the
            # inverse of that. Otherwise rotate to 150 deg (5 * pi / 6)
            target_magnitude = min(abs(angle), 5*math.pi/6)
            rotangle = math.pi - abs(angle) + math.pi - target_magnitude
            if angle2 < 0:
                rotangle *= -1
            hinge = hinge2 - hinge1
            new_coords -= hinge1
            affected = set()
            for _, info in data[corner_idx].items():
                affected = affected.union(info["downstream_indices"])
            for i in affected:
                new_coords[i] = np.dot(rotation_matrix(hinge, rotangle), new_coords[i])
            newcorner = new_coords[corner_idx] + hinge1

            # rotate hinge1
            rotangle = dihedral(ocorner, rod1, hinge1, newcorner) - dihedral(ocorner, rod1, hinge1, corner)
            affected = set()
            for _, info in data[hinge1_idx].items():
                affected = affected.union(info["downstream_indices"])
            for i in affected:
                new_coords[i] = np.dot(rotation_matrix(hinge1-rod1, rotangle), new_coords[i])
            new_coords += hinge1 - hinge2

            # rotate hinge2
            rotangle = dihedral(ocorner, rod2, hinge2, newcorner) - dihedral(ocorner, rod2, hinge2, corner)
            affected = set()
            for _, info in data[hinge2_idx].items():
                affected = affected.union(info["downstream_indices"])
            for i in affected:
                new_coords[i] = np.dot(rotation_matrix(hinge2-rod2, rotangle), new_coords[i])
            new_coords += hinge2
                

    return new_coords


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
