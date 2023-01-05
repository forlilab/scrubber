import numpy as np

def norm(v):
    return v / np.sqrt(np.dot(v, v))

def process_ring6(coords, indices, debug=False):
    assert(len(indices) == 6)
    # for all three pairs of opposing edges, compute their angle

    boat_score = 0.0
    flip_score = [0] * len(indices) # higher score means we should flip atom to convert boat to chair
    dot_edges = [0] * len(indices) # how well aligned are the edges that flip each atom

    for i in range(3):
        a = np.array(coords[indices[(i+0) % 6]])
        b = np.array(coords[indices[(i+1) % 6]])
        c = np.array(coords[indices[(i+3) % 6]])
        d = np.array(coords[indices[(i+4) % 6]])
        edge1 = norm(b-a)
        edge2 = norm(c-d)
        dot = np.dot(edge1, edge2)
        dot_edges[(i+2) % 6] = dot
        dot_edges[(i+5) % 6] = dot

        if debug:
            print(
                "edge1: %2d %2d" % (indices[(i+0)%6]+1, indices[(i+1)%6]+1),
                "edge2: %2d %2d" % (indices[(i+3)%6]+1, indices[(i+4)%6]+1),
                "angle: %.2f" % np.degrees(np.arccos(dot)),
                "dot: %.3f" % dot,
            )

        # up / down
        centroid = np.mean([a, b, c, d], axis=0)
        edge_points = [a, b, c, d]
        corner1 = coords[indices[(i+2) % 6]]
        corner2 = coords[indices[(i+5) % 6]]
        # c1 and c2 are the out of plane magnitudes. They can be either
        # positive or negative depending on direction in which the atom
        # projects out of plane.
        c1 = 0.0
        c2 = 0.0
        # we have four atoms that define two edges, and a centroid.
        # let's iterate over the four possible trianges defined by the centroid
        # and two adjacent atoms, compute the normal for each of the triangles,
        # and check whether the corners are up or down (positive or negative c1/c2)
        for j in range(4):
            v1 = norm(edge_points[(j+1) % 4] - edge_points[j % 4])
            v2 = norm(centroid - edge_points[(j+1) % 4])
            normal = np.cross(v1, v2)
            c1 += np.dot(normal, norm(corner1-centroid))
            c2 += np.dot(normal, norm(corner2-centroid))

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
        if debug:
            print("c1: %.3f" % c1, "c2: %.3f" % c2, "product = %.3f" % (c1*c2))
    return boat_score, flip_score, dot_edges
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
        print(idx)
        other_ring_idxs = [i for i in indices if i != idx] # all ring atoms except the current one
        atom_info = {}
        for neigh in mol.GetAtomWithIdx(idx).GetNeighbors():
            info = {
                "nr_neighbors": 0,
                "downstream_indices": [],
                "bites_own_tail": False,
            }
            neigh_idx = neigh.GetIdx()
            if neigh_idx in other_ring_idxs:  
                continue
            if neigh.GetAtomicNum() == 1:
                continue
            print(idx, neigh_idx, neigh.GetAtomicNum())
            info["nr_neighbors"] = len([a for a in neigh.GetNeighbors() if a.GetAtomicNum() > 1])
            visited = [idx, neigh_idx] 
            info["bites_own_tail"] |= bite_own_tail_recursively(mol, neigh_idx, visited, other_ring_idxs)
            visited.pop(0) # remove ring atom (idx)
            info["downstream_indices"].extend(visited)
            atom_info[neigh_idx] = info
        data[idx] = atom_info
    return data
    
        
             
