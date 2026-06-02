from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem

from rdkit.ForceField.rdForceField import ForceField
import math

import numpy as np
import random

from .openff_toolkit import EspalomaMinimizer

def optimize_conformers(mol: Mol, ff: str="mmff94", max_ff_iter: int = 400):
    optimized_energies = []
    
    mol = Chem.Mol(mol) # create copy

    if (ff == "espaloma"):
        ff = EspalomaMinimizer()
        mol, energies = ff.minimize(mol)
        optimized_energies = list(zip(range(len(energies)), energies))

    else:
        for conf_id in range(mol.GetNumConformers()):


            if ff=="mmff94" and AllChem.MMFFHasAllMoleculeParams(mol):
                ff : ForceField = AllChem.MMFFGetMoleculeForceField(mol, 
                                                            AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94'), 
                                                            confId=conf_id)
            elif ff=="mmff94s" and AllChem.MMFFHasAllMoleculeParams(mol):
                ff : ForceField = AllChem.MMFFGetMoleculeForceField(mol, 
                                                            AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94s'), 
                                                            confId=conf_id)
            else:
                ff : ForceField = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        
            success = ff.Minimize(maxIts=max_ff_iter)
            energy = 0.0
            # print(f"Success: {success}")
            # Check if minimization is successful. 
            if success == 0:
                energy = ff.CalcEnergy()
                optimized_energies.append((conf_id, energy))
            else:
                success = ff.Minimize(maxIts=2*max_ff_iter)
                energy = ff.CalcEnergy()
                optimized_energies.append((conf_id, energy))
    
    return mol, optimized_energies

def add_conformers_to_mol(mol: Mol, conf_coords_list):
    mol = Chem.Mol(mol)  # Make a copy to avoid modifying the original mol
    mol.RemoveAllConformers()  # Clear any existing conformers

    for conf_id, coords in enumerate(conf_coords_list):
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
        conf.SetId(conf_id)
        mol.AddConformer(conf, assignId=True)

    return mol

def find_best_conformer(mol: Mol, ps, num_confs=3, num_etkdg_attempts=1, max_ff_iter=400, ff="mmff94s"):
    """
    Generate multiple conformers with ETKDG and select the one 
    with the lowest energy
    """

    attempts = num_etkdg_attempts

    cids = rdDistGeom.EmbedMultipleConfs(mol, num_confs, ps)

    if mol.GetNumConformers() < 1 and attempts > 1:
        print(f"First ETKDG attempt failed. Will try again until {num_etkdg_attempts} attempts")

    while (mol.GetNumConformers() < 1 and attempts > 1):
        print("attempt: ", num_etkdg_attempts - attempts + 2)
        ps.randomSeed = random.randint(1,100)
        cids = rdDistGeom.EmbedMultipleConfs(mol, num_confs, ps)
        attempts -= 1 

    # if it still fails... 
    if mol.GetNumConformers() < 1:
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unnamed"
        raise ValueError(f"\nETKDG conformer generation failed for molecule: {name} \n"+ 
                         f"Your molecule may be too large or too weird for ETDKG \n"+
                         f"Consider rerunning with higher --num_etkdg_attempts value. \n"+
                         f"If all else fails, you can use --use_random_coords for slower \n"+
                         "but more robust embedding. \n")

    mol, energies = optimize_conformers(mol, ff, max_ff_iter)

    if not energies:
        raise ValueError("No conformers could be optimized during initial generation.")
    
    best_conf_id, best_energy = min(energies, key=lambda x: x[1])


    # Create a new molecule with only the best conformer
    best_mol = Chem.Mol(mol)
    conf = mol.GetConformer(best_conf_id)

    best_mol.RemoveAllConformers()

    best_mol.AddConformer(conf, assignId=True)

    return best_mol, [conf.GetId() for conf in mol.GetConformers()]

#debug
def write_conformers_to_sdf(mol, filename="test.sdf"):
    writer = Chem.SDWriter(filename)
    
    for conf_id in range(mol.GetNumConformers()):
        mol.SetProp("_Name", f"Conformer {conf_id}")  # Optional: Label conformers
        writer.write(mol, confId=conf_id)
    
    writer.close()
    print(f"All conformers written to {filename}")


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
